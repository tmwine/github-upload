function MultAlpha_Iso_Fam(fname)

%{

	function to determine any isomorphisms between block rows of a given Si matrix set

	inputs:
	fname--file name of file containing data_header, Si_mx_cell, and graph_array variables

	outputs:
	br_iso_families--a "nested" cell array; the outer array is a columnar cell vector, with index corresponding to the (increasing-ordered) values in the set of unique values of block row sizes (eg if the block row matrix dimension vector, br_dims, for the data set is [21,16,22,16,21,22,22,16], then row one corresponds to block row matrix dim 16, row two to dim 21, row three to dim 22); the inner cell arrays, referenced by any row in the outer array, contain the isomorphic block row families (a cell array of such grouping vectors), indexed by their (block rows') position in the block row size (br_dims) vector in the data_header data structure
	the option is given to save br_iso_families at the function's end

%}



default_file = ""; % optional; insert name here for default save file


%***
% read in file data
%***

try
	load(fname,"data_header","Si_mx_cell","graph_array");
catch
	printf("problem with input filename in Iso_Fam\n");
	return
end_try_catch

k_sym = data_header.k_sym;
br_dim_vec = data_header.br_dims; % an ordered list of the matrix block row sizes



%***
% main loop, over unique block row matrix sizes
%***

br_sizes = unique(br_dim_vec); % unique also sorts in increasing order
isom_found = false;
tot_fam_found = 0;

for ll=1:length(br_sizes)

	br_ix = find(br_dim_vec==br_sizes(ll));
	printf("block matrix size %i\n",br_sizes(ll));
	num_fams = 0;


	%***
	% find all isomorphism families for all block rows with matrix size br_sizes(ll)
	%***

	while !isempty(br_ix)
		% pick off the first element of br_ix vector, of all block row indices of block rows of size br_sizes(ll), and process the correspondingly-indexed matrix block row; this (block row corresponding to this index) will either match one of the block rows in the existing families, or not; if it does, add it to that family; if not, make a new family

		graph_base = graph_array(br_ix(1));
		mx_filed = false;
		for j1=1:num_fams % see if matrix at block row br_ix(1) fits any existing families
			graph_test = graph_array(sz_fams{j1}(1)); % sz_fams{j1} is a vector, with possibly more than one element in the family; we may as well pick the first, since all within any family are isomorphic by def'n
			isom_found = Iso_Perm_Layer(graph_base, graph_test, data_header.k_sym, data_header.di_vals);
			if isom_found
				sz_fams{j1} = [sz_fams{j1} br_ix(1)];
				mx_filed = true;
				break; % break out of this for loop
			endif
		endfor

		if !mx_filed % this block row isn't isomorphic to any existing families
			num_fams++; % how many families found at this block row matrix dimension (this "ll")
			tot_fam_found++; % how many families found overall
			sz_fams{num_fams} = br_ix(1); % this index doesn't match any existing families, so create a new family
		endif
		br_ix(1) = [];

	endwhile

	br_iso_families{ll} = sz_fams;
	clear sz_fams;

endfor


%***
% option to save, then exit
%***

fflush (stdout);
str_yn = input("Save the results (to an Octave .m file) (y/n)? ","s");
if str_yn!="Y" && str_yn!="y"
	return;
endif

if !strcmp(default_file,"")
	printf("Enter filename (return for default file, %s): ",default_file);
	fname = input("","s");
	if strcmp(fname,"")
		fname = default_file;
	endif
else
	fname = input("Enter filename: ","s");
endif
save(fname,'data_header','br_iso_families');




endfunction







