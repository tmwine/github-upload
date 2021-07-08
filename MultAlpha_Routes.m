function [Si_mx_cell,graph_array] = MultAlpha_Routes(mx_dim,sym_cells,di_vals)


%{

	routine that takes the set of k_sym post-pared binary matrices and finds all closed routes on symbol shifts, separating them into individual block rows, encoded as binary matrices

	inputs:
	mx_dim--dimension of binary matrices to find routes within
	sym_cells--a cell row vector array, each of k_sym entries a 2-column matrix representing, in sparse form, the location of "1"s in the binary shift matrices (of form [r1,c1;r2,c2;...]); the sum of all k_sym of these matrices will have no rows of all zeros and no columns of all zeros
	di_vals--row vector of k_sym entries, representing symbol shift amounts

	outputs:
	Si_mx_cell--a cell array of size [nbrows,k_sym], where nbrows is the number of block rows; the i,j entry contains a 2 column matrix, in format [r1,c1;r2,c2;...], representing the i-th block row's, j-th symbol's binary regularity template matrix
	graph_array--an array of graph structures, nbrows in number, encoding the graphs and supplementary data corresponding to each block row in Si_mx_cell

%}


%***
% setup
%***

% create tot_mix--each row of form [sym#, r, c]
k_sym = length(sym_cells); % number of symbols
tot_mix = [];
tot_mix_nr = 0; % to count total number of rows in tot_mix
for j=1:k_sym
	[r1,c1]=size(sym_cells{j});
	tot_mix = [tot_mix; j*ones(r1,1),sym_cells{j}]; % this will amount to the indicies of "1" values in the sum of all the input symbol shift matrices; rows look like [sym#, r, c], with all symbol #s (1-k_sym) covered
	tot_mix_nr += r1;
endfor

% the chain goes in order col # -> row #
init_col_val = 1; % the initial column value; this amounts to the initial coordinate vector fed into the shift cycle, e_{init_col_val}; we're starting on the 1st column of the total matrix the coordinates in tot_mix represent
cols_avail = [1:mx_dim]; % this is the total shift matrix size; it should be equal to the largest r or c index in the sym_cells paired coords
routes_cell_1 = cell(1,mx_dim); % collects basic [r,c] values along a given route / graph / block row; after further processing, this cell array becomes the central output variable Si_mx_cells; routes should cover all rows (or columns) of the (mx_dim x mx_dim) shift matrices; nb: the cell likely never needs to be this large (and it technically doesn't even need declaring here)
routes_cell_2 = cell(1,mx_dim); % collects the nodes_checked_mx matrices, one from each graph (block row); further on, it is used to create the cyl_graph_matrix fields in output graph_array, the array of graph structures (one for each block row)
num_br_found = 0; % running total of how many block rows (graphs) have been discovered
total_edges_discovered = 0; % counter, for debugging / code error checking



%***
% find and record all routes, based on entry point init_col_val
%***

% this central loop finds all "routes" (or graphs) through the k_sym sym_cell matrices (passed from the paring step); each route (or graph) corresponds to a block row
% the central variables are, new_nodes_last, nodes_found_mx, nodes_checked_mx, and new_nodes_now
	% new_nodes_last are the nodes to check from the previous pass through the loop (or from initialization)
	% nodes_found_mx contains all nodes turned up from following exit edges from all nodes in new_nodes_last
	% nodes_checked_mx contains all nodes that have had their exit edges checked (followed/explored) from previous passes through the loop
	% new_nodes_now contains all nodes in nodes_found_mx that aren't already in nodes_checked_mx
% the main logic through the inner while loop ("while !isempty(new_nodes_last)") is as follows:
	% exploring all nodes in new_nodes_last produces nodes_found_mx
	% the intersection of nodes in nodes_found_mx with those in nodes_checked_mx are the "loopback" nodes (the relevant edges "loop back" to nodes already explored)
	% the remainder of nodes are assigned to new_nodes_now (those not already checked)
	% at the end of the loop, new_nodes_now are assigned to new_nodes_last
	% the loop is repeated, until there are no more nodes to connect to from new_nodes_last
% there are some other things going on in the loop: a total time (tot_time) assigned to each pass through the loop, and assigning indexes [1:n] to the "n" nodes discovered at a given time step, and details on the format of the nodes_found_mx matrix, but the central logic is the backbone; the combination of route finding and graph construction was done all in the same loop because of speed concerns (the code would be clearer if they were separated)

while !isempty(cols_avail) % loops till all permutation routes / graphs found

	% this starts at some r,c "column" value (namely init_col_val) from tot_mix, ie a target value for the "c" of the r,c pairs; it then branches outward, looking for all connections "from" that "c" to all associated "r"'s; the r's then become c's, and it repeats; in this context, the "c" values of tot_mix's rows, [symbol#, r, c], are regarded as node numbers--what gets stored in eg new_nodes_last, and kept as indices in the nodes_checked_mx matrix

	new_nodes_last = [init_col_val,1]; % a central variable, holding the numbers of nodes whose exit edges will be explored this pass through the loop; its rows are, [tot_mix node #, node # at this time step]
	nodes_checked_mx = []; % a central variable; for each graph, this consolidates the node and edge information
	tot_time = 0; % what "time step" we're on in the progression through this graph's node links
	ix_to_remove = []; % keeps track of what rows to remove in tot_mix after the given graph / block row is fully explored
	node_nums_checked = []; % what nodes to remove from cols_avail after the given graph is fully explored (echoes some of the info in nodes_found_mx)
	rte_coor = []; % for construction of the Si matrices; records all the [r,c] pairs discovered while searching for all nodes on a given graph / block row (following this "route" through the [r,c] pairs in tot_mix)
	while !isempty(new_nodes_last) % focus on one total route / graph
		new_tot_ixs = []; % accumulator of indices of tot_mix's rows, corresponding to nodes as they're discovered; used for rte_coor and ix_to_remove near the end of the loop
		tot_time++;
		node_num = 0; % node number at a given time step
		nodes_found_mx = []; % accumulator for node details this time step; this will be checked for loopbacks at the end of this time step, then adjusted for any loopbacks, then appended to nodes_checked_mx
		for j=1:length(new_nodes_last(:,1))

				ix_rows = find(tot_mix(:,3)==new_nodes_last(j,1)); % column vector; of total_mix's rows, [sym#, r, c], this finds any rows in tot_mix (denoted by ix_rows) with a matching column #--in the symbol: col#->row# edge dynamic, this indexes all edges leaving the node enumerated here by the column# denoted in new_nodes_last(j); this should always be non-empty--by the zero-row/-column paring step, every possible should be "covered" in tot_mix; the only time this might fail, is if a previous pass through the route finding loop did not find all connected nodes, and passed ~orphans forward--then we can get "uncovered" nodes
				if isempty(ix_rows)
					"problem with node finding in routes"
					tot_time
					ix_rows
					new_nodes_last
					return
				endif

				new_tot_ixs = [new_tot_ixs; ix_rows]; % for rte_coor

				tmp_len = length(ix_rows);
				nodes_found_mx = [nodes_found_mx; tot_time*ones(tmp_len,1), new_nodes_last(j,2)*ones(tmp_len,1), tot_mix(ix_rows,1), tot_mix(ix_rows,2), new_nodes_last(j,1)*ones(tmp_len,1), (tot_time+1)*ones(tmp_len,1), zeros(tmp_len,1)]; % a central variable; this soon gets split between nodes already in nodes_checked_mx, and those not (which go on to go in new_nodes_now, and then new_nodes_last, for the next pass through the loop)

				total_edges_discovered += length(ix_rows); % for debugging / error checking; this reflects the number of rows "discovered" in tot_mix so far, as we explore all possible graphs; rows in tot_mix should be expressible as a disjoint union of sets, indexed by the values in its third column (the "c" value of the tot_mix rows, [sym#, r, c])--this ensures that's the case (as we loop through all possible "c" values, ie node #s)

		endfor % end loop through new_nodes_last

		node_nums_checked = [node_nums_checked; new_nodes_last(:,1)];

		if !isempty(nodes_checked_mx) % the following will execute if we've already been through this nodes-checking loop (to "prime" nodes_checked_mx); this essentially does the set cutting operation, "cutting" nodes_found_mx's nodes into those that are in nodes_checked_mx already (so they've been covered), and those not in nodes_checked_mx (so they still need to be explored)
			[bool,match_ix] = ismember(nodes_found_mx(:,4),nodes_checked_mx(:,5)); % ~heavy use of ismember here--this will find all instances of values of the 2nd argument vector, in the first; it returns a boolean vector, corresponding to positions in first argument vector that matched a value in 2nd arg. vec., and an index vector, with the index of the matching value in 2nd arg vec at any spots in 1st arg vec that matched; if there are repeated values in the 2nd argument vector, and an element in the 1st arg vec matches those values, the index vector returned will pick just one instance

			% are there loopbacks? ie what nodes in nodes_found_mx are already in nodes_checked_mx
			pos_ix = find(match_ix>0); % which row indices in nodes_found_mx have a matching exit node to those already covered (those in the nodes-covered column of nodes_checked_mx)
			if !isempty(pos_ix) % there are loopbacks
				tmp_1 = nodes_checked_mx(match_ix(pos_ix),:); % this will save rows from nodes_checked_mx where there were matches of its tot_mix node # to the node number in nodes_found_mx
				nodes_found_mx(pos_ix,6:7) = tmp_1(:,1:2); % takes care of loopbacks
			endif
		else % ie 1st pass through the node checking loop (when nodes_checked_mx hasn't yet been started)
			bool = zeros(length(nodes_found_mx(:,4)),1); % every node ref'd by nodes_found_mx(:,4) is a new exit node (this must be the first time step, so there can be no loopbacks)
		endif

		% have new nodes been discovered (ie not all exit nodes were loopbacks)?
		tmp_a = find(bool==0); % row indices in nodes_found_mx of all exit nodes that haven't been covered
		if !isempty(tmp_a) % new nodes have been discovered

			tmp_b = nodes_found_mx(:,4)(tmp_a); % all row indices of nodes_found_mx that have a new exit node
			[new_nodes_now,aa,bb] = unique(tmp_b); % another central variable; new_nodes_now can be thought of in set terms as all nodes in nodes_found that are not in nodes_checked_mx (nodes_found \ nodes_checked_mx); the "unique" removes any redundancies in nodes we've been led to from looping through new_nodes_last--multiple nodes from new_nodes_last could have pointed to the same node # (have the same "r" value, in their respective rows of tot_mix); bb has same dimensions as tmp_b

			nodes_found_mx(tmp_a,7) = bb;
			new_nodes_last = [new_nodes_now,[1:length(new_nodes_now)]']; 

		else % there are no new nodes--all are loopbacks
			new_nodes_last = [];

		endif

		rte_coor = [rte_coor; tot_mix(new_tot_ixs,:)]; % pulling rows from tot_mix corresponding to those indexed in first column of new_nodes_now--so rows will have format [symbol#, r, c]
		ix_to_remove = [ix_to_remove; new_tot_ixs]; % rows to remove from tot_mix, when we're done

		nodes_checked_mx = [nodes_checked_mx; nodes_found_mx];

	endwhile % all possible nodes in this particular graph have been explored

	num_br_found++; % once we've made it to here, what we've just found constitutes another block row

	% rte_coor at this point will contain rows of type [sym#, r, c], all the rows coding the given closed, linked permutation loops for this graph / route

	routes_cell_1{num_br_found} = sortrows(rte_coor); % this will produce [r,c] matrices in Si_mx_cell that are lexicographically ordered; this sorting though adds (processing) time to the routine

	% for binary edge matrices (graph construction)
	routes_cell_2{num_br_found} = nodes_checked_mx;

	% remove the just-found routes from play
	tot_mix(ix_to_remove,:) = [];
	cols_avail = setdiff(cols_avail,node_nums_checked); % this removes all values in node_nums_checked from the available potential node list (cols_avail); note that cols_avail has to be unique for this to work (setdiff can't remove multiple instances of a given value)

	% reset starting column number
	init_col_val=min(cols_avail);

endwhile % end looking for all possible graphs (this forms the block rows)



%***
% post route-finding loop checks and setup
%***

nbrows = num_br_found; % the total value of counter num_br_found, after the loop is done, equals the number of block rows

% debug / error check; this should never trigger--if it does, there's probably a problem in the code; it may be that the above node searching routine is somehow redundantly identifying nodes, or at least accidentally duplicating indices when determining ix_rows (ix_rows should, cumulatively, index every row in tot_mix, exactly once)
if !isempty(tot_mix)
	printf("problem in routes; %i rows left in tot_mix\n",size(tot_mix)(1));
endif
if total_edges_discovered != tot_mix_nr % ie the number of edges discovered should exactly equal the number of rows in tot_mix
	printf("problem in routes; %i edges discovered, out of %i possible edges (values should match)\n", total_edges_discovered, tot_mix_nr);
endif

temporal_graph_cell = cell(nbrows,1); % nb: at this point, the time values can be ~any amount; we have not yet performed the "mod tau" ~cyclical wrap-around / consolidation
for jj=1:nbrows
	temporal_graph_cell{jj} = routes_cell_2{jj}; % includes the original node numbers
endfor

Si_mx_cell = cell(nbrows,k_sym); % for the S_i (sparse form); each cell row will have a set of matrices of [r,c] coordinates (of form [r1,c1;r2,c2;...]) corresponding to the appropriate symbol number (cell column)

bin_val_cell = cell(nbrows,1); % for intermediary computations in upcoming function call for binary edge sequence cell array of matrices



%***
% loop to create the S_i block diagonal matrices; create initial binary edge matrices, for graph generation
%***

% this further parses routes_cell_1 [symbol#,r,c] row-triples into separate S_i matrices (sparse form)
tmp_cell_1 = cell(1,k_sym);
for block_row=1:nbrows

	% consolidate (r,c) pairs by symbol # (ie group into a given matrix)
	for sym_num=1:k_sym
		rcre = find(routes_cell_1{block_row}(:,1)==sym_num);
		tmp_cell_1{sym_num} = routes_cell_1{block_row}(rcre,2:3);
	endfor

	% re-index all matrices on a given cell array row (ie block row), "compressing" the [r,c] values
	all_sc = [];
	for sym_num=1:k_sym
		all_sc = [all_sc; tmp_cell_1{sym_num}];
	endfor
	uc = unique(all_sc); % will be in sorted order; "unique" here is going through all rows and both columns of the [r1,c1; r2,c2; ...] matrix of all_sc
	for sym_num=1:k_sym
		for j1=1:length(uc)
			tmp_cell_1{sym_num}(tmp_cell_1{sym_num}==uc(j1)) = j1;
		endfor
		Si_mx_cell{block_row,sym_num} = tmp_cell_1{sym_num};
	endfor

	% this section is for graph generation purposes, compressing the original tot_mix node numbering in columns 4 and 5 of the respective (block row's) temporal graph matrix, and doing preliminary work for binary degree sequence matrices (here, we assign the sum of the binary edge totals to the resp. (orig) node numbers, for both the "entry" side of the nodes, and the "exit" side
	tmp_mx_1 = temporal_graph_cell{block_row};
	tmp_mx_2 = tmp_mx_1(:,4:5); % extract the two columns with the tot_mix node numbering in them
	pow_sym_vec = 2.^(tmp_mx_1(:,3)-1); % prepares symbol numbers for adding--1->2^0=1; 2->2^1=2; 3->4; 4->8; etc.
	for j1=1:length(uc)
		sym_ix = (tmp_mx_2==uc(j1));
		tmp_mx_2(sym_ix) = j1;
		% these next two lines do the ~edge data (entry and exit) to node bindings
		bin_in(j1) = sum(pow_sym_vec(sym_ix(:,1))); % ie this is of column 4 of temporal_graph_cell--the nodes the resp edges lead in to; this amounts to "counting" (coding as a binary) every edge entering node #j1 (orig, compressed node numbering)
		bin_out(j1) = sum(pow_sym_vec(sym_ix(:,2))); % col 5 of temporal_graph_cell--nodes that the resp edges leave out from; this amounts to "counting" (coding as a binary) every edge leaving node #j1 (orig, compressed node numbering)
	endfor
	tmp_mx_1(:,4:5) = tmp_mx_2; % puts compressed orig node numbers back into tmp_mx_1
	temporal_graph_cell{block_row} = tmp_mx_1;
	bin_val_cell{block_row} = [bin_in;bin_out];

endfor



%***
% for block row graph generation: perform mod tau on the matrices in temporal_graph_cell (through function call); complete binary degree sequences (through function call)
%***

% perform the mod min_orbit on temporal_graph_cell, producing the cylindrical graphs; and finish the binary degree sequence matrices
graph_array = MultAlpha_Graph_Gen(k_sym, di_vals, nbrows, temporal_graph_cell, bin_val_cell);




endfunction







 
