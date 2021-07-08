function out_vec = MultAlpha_Pare(mx_dim,sym_cells)

%{

	routine for removing columns of zeros and rows of zeros from the sum of k_sym binary matrices

	inputs:
	mx_dim--the initial dimension of the binary matrices; it is equal to the total number of unit cells in the overall simplex (of height hk)
	sym_cells--a cell row vector array, each of k_sym entries containing a matrix of unit cell index pairs, corresponding to mappings (of unit cells to other unit cells) under the respective symbol shift (of form [pre-shift index 1,post-shift index 1; pre-shift index 2, post-shift index 2;...])--these pairs may be thought of as "1" entries in a binary matrix

	outputs:
	curr_md--the matrix dimension after the zero-row and zero-column paring step; set to -1 if there is nothing left after paring
	sym_cells_out--a cell row vector array, each of k_sym entries a 2-column matrix representing, in sparse form, the location of "1"s in the binary shift matrices (of form [r1,c1;r2,c2;...]); the sum of all k_sym of these matrices will have no rows of all zeros and no columns of all zeros

%}

%***
% setup
%***

k_val = length(sym_cells);
if k_val < 2
	"problem with input sym_cells"
	return;
endif

% create the total coordinate array (sparse form)
tot_ixs = [];
for j=1:k_val
	[r1,c1] = size(sym_cells{j});
	tot_ixs = [tot_ixs; [sym_cells{j} j*ones(r1,1)]]; % tot_ixs will have [r c sym#] in its rows; sym_cells{j} has # rows, 2 cols, each row a r,c matrix coordinate
endfor



%***
% main paring loop
%***

% this loop makes heavy use of Octave's set manipulation routines: setdiff, ismember, and unique
col_rem = true;
row_rem = true;
nr_ti = size(tot_ixs,1); % number of rows in tot_ixs
curr_md = mx_dim; % current matrix dimension; initially it's the total number of possible unit cells in the entire simplex

while (col_rem || row_rem) && !isempty(tot_ixs)

	ia1 = [];
	ia2 = [];

	zer_rows = setdiff([1:curr_md],tot_ixs(:,1)); % (row vector) gives indices of all-zero rows; careful w/ setdiff--it will not flag repeated values (since [1:curr_md] has no repeats, and we don't care if tot_ixs(:,1) has repeats, it's OK); will return elements in [1:curr_md] (ie row numbers) that are not in tot_ixs(:,1) (ie sparse matrix notation list of rows that have 1's in them)
	if isempty(zer_rows)
		row_rem = false;
	endif
	zer_cols = setdiff([1:curr_md],tot_ixs(:,2)); % (row vector) gives indices of all-zero columns; will return elements in [1:curr_md] (ie col numbers) that are not in tot_ixs(:,2) (ie sparse matrix notation list of cols that have 1's in them)
	if isempty(zer_cols)
		col_rem = false;
	endif

	if row_rem % find indices of all columns corresponding to zer_rows
		tmp_vec_a = [1:nr_ti];
		ia1 = tmp_vec_a(ismember(tot_ixs(:,2),zer_rows)); % (row vector) gives indices in tot_ixs column 2 of values that match any/all values in zer_rows; ie look for which columns match to the zero rows; ia1 will contain the indices of elements in tot_ixs(:,2) (c values of (r,c) pairs) that equal values in zer_rows
	endif
	if col_rem % find indices of all rows corresponding to zer_cols
		tmp_vec_b = [1:nr_ti];
		ia2 = tmp_vec_b(ismember(tot_ixs(:,1),zer_cols)); % (row vector) ia2 will contain the indices of elements in tot_ixs(:,1) (r values of (r,c) pairs) that equal values in zer_cols
	endif

	if row_rem || col_rem
		tot_ixs([ia1';ia2'],:) = [];
		% "compress" coordinates; note this "automatically" removes zero-rows and zero-columns (we compress row/column indices "around" the missing row/col indices, ie those that represent zero-rows and zero-cols)
		nr_ti = size(tot_ixs,1); % refresh row count of tot_ixs, after row deletions
		[dum,dum,tmp_ix] = unique(tot_ixs(:,1:2)); % tmp_ix will be a column vector
		tot_ixs(:,1:2) = reshape(tmp_ix,nr_ti,2); % restore to 2 column matrix
        if !isempty(tot_ixs)
    		curr_md = max(max(tot_ixs(:,1:2))); % update (post-shrinking) matrix dimension
        else
            curr_md = -1; % the matrix set has been pared to nothing
        endif
	endif

endwhile



%***
% recast output as matrices, and return it
%***

sym_cells_out = cell(1,k_val);
for jj=1:k_val
	sym_cells_out{jj} = tot_ixs(tot_ixs(:,3)==jj,1:2);
endfor

out_vec = [curr_md,sym_cells_out]; % this will be converted as a function return value into a cell array of size k_val+1, the "md" value prefixed in its own cell to the cells of sym_cells



endfunction
