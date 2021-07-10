function MultAlpha_Graph_Disp(Si_file,block_row,outfile)

%{

	graph array reader, to create `wraparound' cylinder graph dot file

	inputs:
	Si_file is the filename for the Si_mx_cell file containing the relevant graph_array structure
	block_row is what block row to create the graph of
    outfile is the name for the file to be output via this process (it might create a ~similar-named ~aux file too, or maybe double-writes to outfile)

	output:
	the output file needs to be compiled with graphviz's Kneato, to create a viewable PostScript image; to compile, from the command line:
    	$dot -Kneato -Gsplines=true -Gsep=.3 -Tps outfile -o [ps file out].ps (the -Gsep option puts a "bubble" around the nodes as drawn, to put whitespace between spline and node)

%}



%***
% create initial graph, and save to a file
%***

Y_TOL = 10; % how much drift to allow average y coord of node columns
INPUT_SCALE = 72; % for converting from points to inches positions; if nodes in neato output are too large / crowded, then increase this value (baseline is 72)

load(Si_file,"data_header","graph_array");

for jj=1:data_header.k_sym
	str_lab{jj} = char('a'+(jj-1)); % for edge labels
endfor

tm_nd_mx = graph_array(block_row).cyl_graph_matrix;

% create a tmp file, "outfile" w/ a "_aux" suffix
tmp_str = strsplit(outfile,".");
outfile_aux = [tmp_str{1}, "_aux"];
for jj=2:length(tmp_str)
    outfile_aux = [outfile_aux, ".", tmp_str{jj}];
endfor

fid = fopen(outfile_aux,'w');
fdisp(fid,"digraph {");
fdisp(fid,"rankdir=\"LR\"");

% add the wraparound column nodes
tot_time = tm_nd_mx(end,1);
max_nd_num = max(tm_nd_mx(:,5)); % maximum ID# used for all nodes in original graph
tmp_ix = find(tm_nd_mx(:,1)==tot_time); % row indices of entries in tm_nd_mx at last time step
tmp_org_ID_vec = tm_nd_mx(tmp_ix,4); % all node IDs of nodes that edges leaving nodes at last time step enter--such nodes must be at the first time step (only); may be repeats (different edges entering same node)
[yy,aa,bb] = unique(tmp_org_ID_vec);
zz = max_nd_num + [1:length(yy)];
wrp_nd_pr = -ones(length(yy),2);
for jj=1:length(tmp_ix)
    if (wrp_nd_pr(bb(jj))==-1)
        wrp_nd_pr(bb(jj),1) = tm_nd_mx(tmp_ix(jj),4); % node ID from node in lh column
        wrp_nd_pr(bb(jj),2) = zz(bb(jj)); % matching node ID of new node for new "wraparound" rh column
    endif
    tm_nd_mx(tmp_ix(jj),4) = zz(bb(jj)); % change the values in tm_nd_mx
endfor

for mm=1:2
    for jj=1:size(wrp_nd_pr,1)
        if (mm==2)
            str_app = ", style=\"dashed\"]";
        else
            str_app = "]";
        endif
        str_out = [num2str(wrp_nd_pr(jj,mm)), " [label=\"", num2str(jj), "\"", str_app];
        fdisp(fid,str_out);
    endfor
endfor
fdisp(fid,"node[label=\"\"];");
tot_nodes = length(unique([tm_nd_mx(:,4);tm_nd_mx(:,5)]));

for jj=1:size(tm_nd_mx,1)
	fdisp(fid,[num2str(tm_nd_mx(jj,5)),"->",num2str(tm_nd_mx(jj,4)),"[label=\"",str_lab{tm_nd_mx(jj,3)},"\"]"]);
endfor

nat_cell = cell(1,tot_time); % for cell of vectors, one for each time step
rank_cell = cell(1,tot_time+1); % to save which nodes go in columns 1,2,...,tot_time+1
for jj=1:tot_time
	tmp_ix = find(tm_nd_mx(:,1)==jj); % 1st col of tm_nd_mx is the time step #
	rank_vec = unique(tm_nd_mx(tmp_ix,5)); % will be a column vector
    rank_cell{jj} = rank_vec; % node ID numbers in column jj
	rank_str = ["{rank=same; ", num2str(rank_vec(1))];
	for kk=2:length(rank_vec)
		rank_str = [rank_str, ", ", num2str(rank_vec(kk))];
	endfor
	rank_str = [rank_str, "}"];
	fdisp(fid,rank_str);
	%tmp_str = num2str(rank_vec);
	%tmp_str_2 = strrep(tmp_str,"  "," ");
	%fdisp(fid,["{rank=same; ", nat_cell{jj},"}"]);
endfor

% group all nodes in last (added) wraparound column into a rank
rank_str = ["{rank=same; ", num2str(wrp_nd_pr(1,2))];
for jj=2:size(wrp_nd_pr,1)
    rank_str = [rank_str, ", ", num2str(wrp_nd_pr(jj,2))];
endfor
rank_str = [rank_str, "}"];
rank_cell{tot_time+1} = wrp_nd_pr(:,2);
fdisp(fid,rank_str);

fdisp(fid,"}");

fclose(fid);



%***
% run "dot -Tdot" on the initial graph file, via system call
%***

str_in = ["dot -Tdot ", outfile_aux, " -o ", outfile];
[status, output] = system(str_in);
if (status != 0)
    "problem with calling dot; check to be sure graphviz / dot is installed"
    return
endif



%***
% scrape the augmented text .dot file for the node position data, using graphviz's gvpr command
%***
str_in = ["gvpr 'N{print ($.name, \",\", $.pos)}' " outfile];
[status, output] = system(str_in);
if (status != 0)
    "problem with calling gvpr; check to be sure graphviz / dot is installed"
    return
endif

out_parse = strsplit(output,'\n');
for kk=1:length(out_parse)
	if !isempty(out_parse{kk})
		tmp_parse = strsplit(out_parse{kk},',');
		if length(tmp_parse)!=3
			"problem with scraping node positions from intermediary 'outfile'"
			return
		endif
		for jj=1:3
			pos_mx(kk,jj)=str2num(tmp_parse{jj}); % pos_mx has rows of format [node #, x_pos, y_pos]; x,y coords have lower left corner as origin
		endfor
	endif
endfor

if (size(pos_mx,1)!=tot_nodes)
    "something wrong with -Tdot and annotated graph file; not enough node position data"
    return
endif



%***
% modify the position data; adjust graph to be more or less horizontal; line up the right column correctly with the left
%***

% organize node column position information; augment rank_cell vectors with position data (create matrices)
for jj=1:tot_time+1
    tmp_vec = rank_cell{jj};
    tmp_mx = [];
    for kk=1:length(tmp_vec) % go through nodes in column jj; rank_cell has column vectors
        tmp_ix = find(pos_mx(:,1)==tmp_vec(kk));
        if (length(tmp_ix) != 1)
            "problem with pos_mx data; node not found, or multiple nodes with same ID"
            return
        endif
        tmp_mx(kk,:) = pos_mx(tmp_ix,:); % node #, x, y
    endfor
    rank_cell{jj} = tmp_mx;
endfor

% adjust for vertical drift in graph node columns (should be fairly level)
for jj=1:length(rank_cell)
    mean_vec(jj) = mean(rank_cell{jj}(:,3)); % avg y coord for node column jj
endfor
tot_avg_ht = mean(mean_vec);
min_y_coord = 0;
for jj=1:length(rank_cell)
    ht_dff = mean_vec(jj)-tot_avg_ht;
    if (abs(ht_dff)>Y_TOL)
        rank_cell{jj}(:,3) -= ht_dff-Y_TOL;
    endif
    tmp_min = min(rank_cell{jj}(:,3));
    if (tmp_min < min_y_coord)
        tmp_min = min_y_coord;
    endif
endfor
if (min_y_coord<0) % if any of the y coords have shifted < 0, shift them all so min is zero
    for jj=1:length(rank_cell)
        rank_cell{jj}(:,3) += abs(min_y_coord);
    endfor
endif

% set positions of rightmost "wraparound" column to match those in lh column (node pair matchings are stored in wrp_nd_pr)
for jj=1:size(wrp_nd_pr,1)
    tmp_ix_1 = find(rank_cell{1}(:,1)==wrp_nd_pr(jj,1)); % row index in lh column's listing in rank cell belonging to jjth node in wrp_nd_pr's 1st col of node pairs (of node IDs of 1st col of graph)
    tmp_ix_2 = find(rank_cell{end}(:,1)==wrp_nd_pr(jj,2)); % row index in rh column's listing in rank_cell belonging to jjth node listed in wrp_nd_pr's 2nd column of node pairs (of node IDs of last col of graph)
    if (length(tmp_ix_1)!=1 || length(tmp_ix_2)!=1)
        "problems in rank_cell listings; node missing, or multiple nodes with same ID"
        return
    endif
    rank_cell{end}(tmp_ix_2,3) = rank_cell{1}(tmp_ix_1,3); % match the y coords for lh/rh col node pair with resp IDs wrp_nd_pr(1), wrp_nd_pr(2)
endfor



%***
% create a new .dot graph text file for use with Kneato
%***
fid = fopen(outfile,'w'); % this will overwrite the appended file from -Tdot--once positions have been scraped it's no longer needed
fdisp(fid,"digraph {");
fdisp(fid,["inputscale=",num2str(INPUT_SCALE)]);
fdisp(fid,"rankdir=\"LR\"");

% force node positions
for jj=1:length(rank_cell)
    tmp_mx = rank_cell{jj}; % rows with [node #, x, y]
    if (jj==1 || jj==length(rank_cell))
        [dum,tmp_ord] = sort(tmp_mx(:,3),"descend");
        ord_num(tmp_ord) = [1:length(tmp_ord)]; % invert to obtain correct node ranking order, by y coordinate
    endif
    if (jj==2)
        fdisp(fid,"node[label=\"\"];"); % don't want to number middle nodes
    endif
    for kk=1:size(tmp_mx,1)
        if (jj==1)
            tmp_str = [num2str(tmp_mx(kk,1)), "[pos=\"", num2str(tmp_mx(kk,2)), ",", num2str(tmp_mx(kk,3)), "!\", label=\"", num2str(ord_num(kk)), "\"];"];
        elseif (jj==length(rank_cell))
            tmp_str = [num2str(tmp_mx(kk,1)), "[pos=\"", num2str(tmp_mx(kk,2)), ",", num2str(tmp_mx(kk,3)), "!\", label=\"", num2str(ord_num(kk)), "\", style=\"dashed\"]"];
        else
            tmp_str = [num2str(tmp_mx(kk,1)), "[pos=\"", num2str(tmp_mx(kk,2)), ",", num2str(tmp_mx(kk,3)), "!\"];"]; % eg: 1[pos="55,162!"];
        endif
        fdisp(fid,tmp_str);
    endfor
endfor
fdisp(fid,"node[label=\"\"];");

% describe edges
for jj=1:size(tm_nd_mx,1)
	fdisp(fid,[num2str(tm_nd_mx(jj,5)),"->",num2str(tm_nd_mx(jj,4)),"[label=\"",str_lab{tm_nd_mx(jj,3)},"\"]"]);
endfor

% describe same-rank nodes (ie node columns)
for jj=1:length(rank_cell)
    rank_vec = rank_cell{jj}(:,1);
	rank_str = ["{rank=same; ", num2str(rank_vec(1))];
	for kk=2:length(rank_vec)
		rank_str = [rank_str, ", ", num2str(rank_vec(kk))];
	endfor
	rank_str = [rank_str, "}"];
	fdisp(fid,rank_str);
endfor

fdisp(fid,"}");

fclose(fid);






endfunction
