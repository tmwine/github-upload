function total_match_found = Iso_Full_Layer(graph_base,graph_test,k_sym,perm_vec,circ_perm_ok,row_lh_start,lh_mx_out_1)

%{

	this function completes the isomorphism checks by applying the received symbol permutation to the "test" binary degree sequence matrix's rows (the edges in/out binary strings), then tries to match the "test" b.d.s. matrix to the "base" b.d.s. matrix (bin_ds_mx_1) via temporal circular permutations; if there's a match it then calculates the left-hand graph of the (symbol-permuted, and circularly rotated) test graph--if this matches the base graph's left-hand graph (lh_mx_out_1), it's found an isomorphism; this function is called once per symbol permutation

	inputs:
	graph_base--the graph structure being tested against
	graph_test--the graph structure whose symbols get permuted, checking at each permutation for isomorphism to the graph represented by graph_base
	k_sym--number of symbols
	perm_vec--the vector of permuted symbol numbers (some permutation of [1,...,k_sym])
	circ_perm_ok--a vector of the temporal offsets, applied to the test graph, where there's any chance of a match; these positions are the only ones it's necessary to examine further--ie for the binary edge in/out row-strings
	row_lh_start--the index (here, the row #, corresponding to rows of bin_edge_mx_base) that was chosen as optimal (of lowest binary edge in/out value multiplicity) of the base graph during the initial S_i graph processing step; the multiplicity can still be >1
	lh_mx_out_1--a coding of the base graph, rendered via left-hand path search, starting at min_multiplicity_node; the final check for isomorphism relies on making an educated guess at the corresponding node in the test graph--if the left-hand graph rendering of the test graph at this node starting point is the same as lh_mx_out_1, then they're isomorphic

	outputs:
	total_match_found--boolean

	nb: this requires / expects both the binary degree sequence matrices to have their paired (timestamp) columns' rows sorted (lexicographical sort)

%}


%***
% setup
%***

total_match_found = false;
perm_match_found = false;
min_orbit = graph_base.min_orbit;
max_node_count = graph_base.max_node_count;
tot_pos = min_orbit*max_node_count; % total number of rows in the node matrices (cyl_graph_matrix matrices, and node_numbering_vec vectors; such rows include, if needed, padding 0's) (used as total shift positions under circular permutations)
bin_edge_mx_base = graph_base.binary_edge_matrix;
bin_edge_mx_test = graph_test.binary_edge_matrix;
test_nn_vec = graph_test.node_numbering_vec;



%***
% change binary edge strings in bin_edge_mx_test matrix to reflect the symbol permutations
%***

% perform permutation and re-sort nodes/rows binary symbol edge assigments within each time step (eg k_sym=3, [in,out] string row is "011100", and under a->b, b->c, c->a symbol permutation (perm_vec=[2,3,1]), row becomes "101010")
if !isequal(perm_vec,[1:k_sym])
	bin_edge_mx_test(:,[perm_vec, k_sym+perm_vec]) = bin_edge_mx_test; % recall the k_sym leftmost cols are for the coding the "in" edges, and the rightmost k_sym cols for coding the "out" edges

	for jj=1:min_orbit % loop through each time step to sort nodes/rows at that time step
		blk_str = bin_edge_mx_test(max_node_count*(jj-1)+(1:max_node_count),:); % the block (consec rows) of strings corresp to this time step (jj)
		[blk_srt,srt_ix] = sortrows(blk_str); % sortrows goes in ascending order
		bin_edge_mx_test(max_node_count*(jj-1)+(1:max_node_count),:) = flipud(blk_srt); % flipud--we want descending order
		tmp_blk_onn = test_nn_vec(max_node_count*(jj-1)+(1:max_node_count));
		test_nn_vec(max_node_count*(jj-1)+(1:max_node_count)) = flipud(tmp_blk_onn(srt_ix)); % the original node numbering assignments need to be permuted just as the sorted binary string rows are at this time step
	endfor

else
	% intentionally left bank; nb: by design (eg routes and graph_gen routines) the incoming binary degree sequence matrices should have their string rows (node binary edge strings) within a given time period already sorted by default (in lexi. descending order); so for the identity permutation, no sorting should be needed
endif



%***
% look for any matches between the test and base binary degree sequence graphs, under circular rotations of the (post-symbol-permuted) test b.d.s. graph
%***

% try to match the binary degree sequence graphs under circular time rotations
perm_match_found = false;
rot_base_bds = bin_edge_mx_base; % the "base" binary degree sequence will remain the same; the "test" (permuted) binary degree sequence will be circularly rotated
circ_match_vec = [];
for jj=1:length(circ_perm_ok) % recall, circ_perm_ok found the temporal offset amounts for "rough" matches between the node count vectors (ie how many nodes at each time step); these positions are the only ones it's necessary to examine further--ie for the binary edge in/out row-strings
	cp_pos = circ_perm_ok(jj); % will be a value in [0,min_orbit-1]
	rot_vec = mod(max_node_count*(min_orbit-cp_pos)+[0:(tot_pos-1)],tot_pos) + 1;
	rot_test_bds = bin_edge_mx_test(rot_vec,:); % "rotates" the rows by max_node_count times cp_pos
	if isequal(rot_base_bds,rot_test_bds)
		perm_match_found = true;
		circ_match_vec = [circ_match_vec, cp_pos]; % store the shift position(s) where a match is found between the binary degree sequence matrices
	endif
endfor

if !perm_match_found
	total_match_found = false;
	% DEBUG
	%printf("no match found at permutation ");
	%printf(num2str(perm_vec));
	%printf("\n");
	return
endif



%***
% check left-hand graph renderings of the test graph, finding the appropriate left-hand rendering start node(s) from the circular rotation match(es), and the corresponding start node(s) from the base graph
%***

% if we've gotten here, there's been a match between the binary degree sequence, in/out edge matrices, under the symbol permutation on the test graph, and under some (possibly trivial) circular time step rotation on the test graph; the last, most definitive step, is to check the left-hand path renderings
% row_lh_start is the index (here, the row #, corresponding to rows of bin_edge_mx_base) that was chosen as optimal (of lowest binary edge in/out value multiplicity) during the initial S_i graph processing step; the multiplicity can still be >1
for jj=1:length(circ_match_vec)


	%***
	% find nodes in the test graph that match the left-hand graph starting node in base graph
	%***

	cp_pos = circ_match_vec(jj);
	row_test = mod((row_lh_start-1)-max_node_count*cp_pos,tot_pos) + 1; % this shifts the circular permutation backwards, to find the corresponding row in the "test" binary degree sequence (bin_edge_mx_test) that lines up with the target row to start from in bin_edge_mx_base
	str_start = bin_edge_mx_test(row_test,:);

	% DEBUG
	if !isequal(str_start,bin_edge_mx_base(row_lh_start,:))
		"error in circular permutations"
		return
	endif

	% DEBUG
	printf("circ perm equal at %i, at symbol permutation  ",cp_pos);
	printf(num2str(perm_vec));

	tmp_ts = ceil(row_test/max_node_count); % determine the time step this row/node is in
	tmp_blk = bin_edge_mx_test(max_node_count*(tmp_ts-1)+(1:max_node_count),:); % block of all nodes (strings) at this time step (including padding zero strings)

	nd_tt = []; % saves the node numbers (original numbering) of nodes to test the left hand graph renderings of
	for mm=1:max_node_count % find all strings at this time step equal to the target string, and save the corresponding original node numbers (of the "test" graph)
		if isequal(tmp_blk(mm,:),str_start)
			nd_tt = [nd_tt,test_nn_vec(max_node_count*(tmp_ts-1)+mm)];
		endif
	endfor

	
	%***
	% with test graph nodes in nd_tt, perform left-hand graph renderings starting at those nodes, and check for a match to the base graph's left-hand rendering
	%***

	for mm=1:length(nd_tt)
		tmp_symcol = graph_test.cyl_graph_matrix(:,3); % column of symbols (1 to k_sym) from the test cyl_graph_matrix
		tmp_newcol = zeros(length(tmp_symcol),1);
		for ll=1:length(tmp_symcol)
			tmp_newcol(ll) = perm_vec(tmp_symcol(ll));
		endfor
		tmp_node_mx = graph_test.cyl_graph_matrix;
		tmp_node_mx(:,3) = tmp_newcol;
		lh_mx_out_2 = Iso_LH_Graph_Chk(k_sym,tmp_node_mx,nd_tt(mm));
		if isequal(lh_mx_out_1(1:k_sym,:),lh_mx_out_2(1:k_sym,:)) % recall, as of latest, the output of the left hand walk function has, on the bottommost (k_sym+1 th) row, the original node numbers--these are not part of the comparison
			printf("; isomorphism found");
			%printf(num2str(perm_vec));
			total_match_found = true;

			% DEBUG
			printf("\n");
			return

		endif
	endfor

	printf("\n");

endfor





endfunction
