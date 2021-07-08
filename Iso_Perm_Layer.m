function isom_found = Iso_Perm_Layer(graph_base, graph_test, k_sym, di_vals)


%{ 

	permutation layer for isomorphism checker
	after some preliminary checks, this cycles through acceptable permutations on symbols (symbols with the same di_vals shift value are permuted, within a "pool"; there may be zero or more such pools in the di_vals shift values vector)
	at each symbol permutation (starting with the identity permutation), it checks for equal symbol (edge) proportions between base and test graphs, and if they match, it calls the full isomorphism checker
	this returns as soon as it's found any isomorphism, at any permutation

	inputs:
	graph_base--the graph structure being tested against
	graph_test--the graph structure whose symbols get permuted, checking at each permutation for isomorphism to the graph represented by graph_base
	k_sym--number of symbols
	di_vals--vector of symbol shift amounts

	output:
	isom_found--boolean; true if any isomorphism has been found

%}

isom_found = false; % this will be true if there is an isomorphism found between the two graphs (after any allowable symbol permutations)
lh_base_graph_created = false;


%***
% preliminary checks; pre-permutation
%***

	% check for time period match
min_orbit = graph_base.min_orbit;
if !isequal(min_orbit,graph_test.min_orbit,graph_base.cyl_graph_matrix(end,1),graph_test.cyl_graph_matrix(end,1))
	% minimum periods don't match (or something is wrong with the last timestamps between the two graphs)
	% DEBUG
	"periods don't match"
	isom_found = false;
	return;
endif

	% check for node count match (matrix dimensions should match; this check may be redundant)
num_nodes = length(unique(graph_base.cyl_graph_matrix(:,5)));
if !isequal(num_nodes,length(unique(graph_test.cyl_graph_matrix(:,5))),sum(graph_base.node_count_vec),sum(graph_test.node_count_vec))
	% graphs don't have the same number of nodes
	% DEBUG
	"different node amounts"
	isom_found = false;
	return;
endif

	% check for consistency between the number of rows in the binary edge sequence matrices and the values of min_orbit and the max node count over all time steps
max_node_count = graph_base.max_node_count;
tmp_nr_1 = size(graph_base.binary_edge_matrix)(1);
tmp_nr_2 = size(graph_test.binary_edge_matrix)(1);
if !(max_node_count==graph_test.max_node_count) || !isequal(tmp_nr_1,tmp_nr_2,min_orbit*max_node_count)
	% binary edge sequence matrices row counts don't add up
	% DEBUG
	"different binary edge sequence parameters"
	isom_found = false;
	return
endif

tot_pos = max_node_count*min_orbit; % total number of rows in the node matrices (binary_edge_matrix matrices, and node_numbering_vec vectors; such rows include, if needed, padding 0's) (used as total shift positions under circular permutations)

	% find any matches under circular permutations of the node count vectors (node_count_vec, which holds the node counts at each time step) between graph_test and graph_base graph structures (these node count values are independent of symbol permutations, and at least one match should be found in order to continue); if there are any matches, save them; if not return without a match
	% the circular permutations are "unit time rightward (or downward) shifts", mod max_node_count*min_orbit, of the involved indices; a shift of zero is the identity permutation
circ_perm_ok = [];
for jj=0:min_orbit-1
	tmp_pix = mod((min_orbit-jj)+[0:min_orbit-1],min_orbit)+1;
	if isequal(graph_test.node_count_vec(tmp_pix),graph_base.node_count_vec)
		circ_perm_ok = [circ_perm_ok,jj];
	endif
endfor
if isempty(circ_perm_ok)
	% no way to match the temporal node counts (these are independent of symbol permutations)
	% DEBUG
	%"temporal node counts won't match"
	isom_found = false;
	return
endif



%***
% initial preparation for symbols permutation loop; determining pools of equal symbol shift values
%***

	% initialize edge count vectors, for more preliminary checks at this layer after each permutation
edge_vec_1 = graph_base.cyl_graph_matrix(:,3);
edge_vec_2 = graph_test.cyl_graph_matrix(:,3);
for jj=1:k_sym
	sc_vec_1(jj) = sum(edge_vec_1==jj); % save "base" graph's edge count vector
endfor

	% find the row in the "base" binary edge matrix corresponding to the preferred (lowest multiplicity of binary edge in/out string codes) node to produce the left hand graph from
row_lh_start = find(graph_base.node_numbering_vec==graph_base.min_multiplicity_node);
if length(row_lh_start)>1 || length(row_lh_start)==0
	% if this triggers, there's a problem with the code; node_numbering_vec should have exactly one occurrence of every node number (along with possible padding 0's)
	"problem with min multiplicity node"
	return
endif

	% parse di_vals into pools of equal values
[udi,dum,u_ix] = unique(di_vals);
for jj=1:length(udi)
	di_pool_size(jj) = sum(di_vals==udi(jj));
endfor
pool_ix = find(di_pool_size>1);
if !isempty(pool_ix)
	num_pools = length(pool_ix);
else
	num_pools = 0;
endif



%***
% if there are no pools (no symbols to permute) just check for identity permutation isomorphism, and return
%***

if num_pools==0 % if there are no permutations to perform, check symbol proportions, then go on to next layer
	% preliminary test for agreement in edge (symbol) proportions

	% DEBUG
	printf("no perms; only identity perm");

	for jj=1:k_sym
		sc_vec_2(jj) = sum(edge_vec_2==jj); % obtain the counts of the test graph's symbol edge vector
	endfor
	if !isequal(sc_vec_1,sc_vec_2)
		% symbol counts don't match between graphs; ie the amounts of edges of each symbol don't agree for this graph pair

		% DEBUG
		printf("; perm symbol counts at identity permutation don't match\n");

		isom_found = false;
		return
	else
		% so far so good; call next (fullest) isomorphism check layer
		if !lh_base_graph_created
			lh_mx_out_1 = Iso_LH_Graph_Chk(k_sym,graph_base.cyl_graph_matrix,...
				graph_base.min_multiplicity_node);
			lh_base_graph_created = true;
		endif
		isom_found = Iso_Full_Layer(graph_base,graph_test,k_sym,...
			perm_vec,circ_perm_ok,row_lh_start,lh_mx_out_1);
		if isom_found
			% DEBUG
			%"match found"
			%perm_vec
			printf("; isomorphism found\n");
			return

		endif
	endif
	printf("\n");
	% if its gotten here, it passed the symbol edge count match, but not the full isomorphism check
	return
endif



%***
% final permutation prep; set up incremental permutation variables
%***

% pool_ix_cell will be a cell array of size equal to the number of nontrivial (size>1) pools in di_vals, the jth vector in the cell array holding indices in di_vals corresponding to the respective (repeated) shift amount, of udi(pool_ix(j))
pool_ix_cell = cell(1,num_pools);
for jj=1:num_pools
	pool_ix_cell{jj} = find(di_vals==udi(pool_ix(jj)));
endfor

% create the cell array for incremental permutation for all pools
perm_inc_cell = cell(num_pools,2);
for j1=1:num_pools
	num_set = pool_ix_cell{j1};
	perm_slots = length(num_set);
	in_mx = [num_set;ones(1,perm_slots)];
	sel_cell = cell(1,perm_slots);
	for j2=1:perm_slots
		sel_cell{j2} = num_set(j2:perm_slots);
	endfor
	perm_inc_cell{j1,1} = in_mx;
	perm_inc_cell{j1,2} = sel_cell;
endfor

% initialize for the main permutation loop
bin_pv = zeros(1,num_pools); % binary pointer vector--tells which pool to increment; eg [0,0,0] indicates increment highest-indexed pool (#3); [0,0,1] indicates increment next-highest-indexed pool (#2); when this pointer vector is "full" ([1,1,1]) it signals all permutations are done
save_inc_cell = perm_inc_cell; % saves, for resets; in the loop, perm_inc_cell is manipulated
inc_flag = false;
all_done = false;



%***
% permutation loop
%***

% the loop starts with the identity permutation, checking for isomorphism, then continues by incrementally permuting the symbols within the one or more pools
% the permutation vector perm_vec reads as symbol mapping [1:k_sym]->perm_vec--eg for a,b,c=1,2,3, perm_vec=[2,3,1] is the mapping a->b, b->c, c->a
% the permuting itself is done incrementally (no tables created beforehand), through the combination of the sel_cell cell array, and the permute_1 function--there are probably better ways to do this though
% nb: once the relevant in_mx 2-row permutation matrix, in the relevant cell of perm_inc_cell has run through all permutations, it requires another call to permute_1 to determine all permutations are done

while !all_done

	% which pool to increment; this increments the highest indexed pool (eg 3 pools, the pool indexed by "3") until it is finished with all permutations, then the next highest once (unless its finished w/ all perms), then the highest again till its finished, etc.
	pti_ix = max(find(bin_pv==0));

	if !inc_flag % if this is a loop pass when inc_flag is true, nothing has changed in the permutation matrices; so don't call any of the isomorphism checkers this pass

		% create the permutation vector, out of the pools
		perm_vec = [1:k_sym]; % symbols (values) at indices matching the pools' (in pool_ix_cell's vectors of indices) get permuted; the others remain fixed
		for jj=1:num_pools
			perm_vec(pool_ix_cell{jj}) = perm_inc_cell{jj,1}(1,:); % this is what the old symbols [1:k_sym] map to, in order--eg old: [1,2,3] and perm_vec=[3,1,2] means 1->3, 2->1, 3->2 for old->new mappings
		endfor

		% preliminary test for agreement in edge (symbol) proportions
		for jj=1:k_sym
			sc_vec_2(jj) = sum(edge_vec_2==perm_vec(jj)); % obtain the counts of the test graph's permuted-symbol edge vector
		endfor
		if !isequal(sc_vec_1,sc_vec_2)
			% symbol counts don't match between graphs; ie the amounts of edges of each symbol don't agree for this permutation
			%DEBUG
			%printf("perm symbol counts don't match, at permutation ");
			%printf(num2str(perm_vec));
			%printf("\n");
		else
			% so far so good; call next (fullest) isomorphism check layer
			if !lh_base_graph_created
				lh_mx_out_1 = Iso_LH_Graph_Chk(k_sym,graph_base.cyl_graph_matrix,...
					graph_base.min_multiplicity_node);
				lh_base_graph_created = true;
			endif
			isom_found = Iso_Full_Layer(graph_base,graph_test,k_sym,...
				perm_vec,circ_perm_ok,row_lh_start,lh_mx_out_1);

			if isom_found
				% DEBUG
				%"match found"
				%perm_vec
				return

			endif
		endif

	endif


	%***
	% update (increment) the pool permutations
	%***

	% call the permutation function, to incrementally permute the symbol-pool indices
	[out_mx,tmp_cell,fin_perm] = permute_1(perm_inc_cell{pti_ix,1},perm_inc_cell{pti_ix,2});

	% if the permutation function returns "done" (i.e. fin_perm==true) increment the binary pointer vector (eg pointer vec at top of this pass through loop was [0,0,0] (necc incrementing the highest indexed pool) we get a "done" (in this case on the highest indexed pool), so we increment to [0,0,1]; if we came in with [0,0,1] (necc incrementing the 2nd highest pool), and got a "done," we'd then have [0,1,1]); and if the binary pointer vector increments to [1,...,1], then leave loop--all perms done
	% if permutation function does not return fin_perm==true and if we incremented last step (flag), then do the reset (all lower pools get perms reset, and pointer vec gets reset)
	if fin_perm
		if pti_ix == 1
			% this means we've just finished the permutations at the lowest indexed (most outer) pool; so we're done with all the permutations in all pools
			all_done = true;
		else
			bin_pv(pti_ix) = 1;
			if !inc_flag
				inc_flag = true;
			endif
		endif
	else
		perm_inc_cell{pti_ix,1} = out_mx;
		perm_inc_cell{pti_ix,2} = tmp_cell;
		if inc_flag % we're here because the pool given by bin_pv pointer vector can be further permuted--if in addition we've just come in from a binary pointer vector increment in the last step, then reset all more inner pools (higher indexed); pti_ix in this case would be the index of the pool just (successfully) further permuted, so all higher indexed pools need to be reset
			for jj=pti_ix+1:num_pools
				perm_inc_cell{jj,1} = save_inc_cell{jj,1};
				perm_inc_cell{jj,2} = save_inc_cell{jj,2};
			endfor
			bin_pv = zeros(1,num_pools);
			inc_flag = false;
		endif
	endif

endwhile




endfunction







function [out_mx,sel_cell,fin_perm] = permute_1(in_mx,sel_cell)

	% checks which counters are full (a full counter at the ith index (i in [1,perm_slots]) has value perm_slots-i+1; reports the first non-full counter, what is the basis for updating the distribution of sets in sel_cell
	fin_perm = false;
	perm_slots = size(in_mx)(2);
	perm_val_vec = [perm_slots:-1:1];
	full_vec = (in_mx(2,:)==perm_val_vec);
	nf_ix = find(full_vec==0);
	if isempty(nf_ix)
		fin_perm = true;
		out_mx = in_mx;
		return
	endif

	% update counters; according to highest (most leftward) full index, reset all lower indices
	ix_to_inc = max(nf_ix); % all places to the right of this are full (complete w/ their permutation cycle); so this index gets the next value, and all index counts to its right are reset (to 1)
	in_mx(2,ix_to_inc)++;
	in_mx(2,ix_to_inc+1:perm_slots) = 1;

	% update sel_cells
	nv = sel_cell{ix_to_inc}(in_mx(2,ix_to_inc)); % this is the counter we just incremented--it now points to the next value in the resp selection set in sel_cell; choosing this value now resets selection sets at all later cells (higher index number) (they can't have this value in them)
	sel_rem = setdiff(sel_cell{ix_to_inc},nv);
	for jj=ix_to_inc+1:perm_slots
		%ss_ix = perm_slots-jj+1;
		sel_cell{jj} = sel_rem((jj:perm_slots)-ix_to_inc);
	endfor

	%sel_cell

	% DEBUG
	%{
	mx_tmp = [];
	for jj=1:perm_slots
		mx_tmp = [mx_tmp; zeros(1,jj-1),sel_cell{jj}];
	endfor
	[mx_tmp; in_mx(2,:)]
	%}

	% "write" sel_cells, according to indexes in in_mx(2,:), to permutation vector to be output (out_mx(1,:))
	for jj=ix_to_inc:perm_slots
		in_mx(1,jj) = sel_cell{jj}(in_mx(2,jj));
	endfor

	out_mx = in_mx;

endfunction



