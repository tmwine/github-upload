function res_o = MultAlpha_Mx_Gen(k_sym, di_vals, hk, delta, p_and_s)

%{

	function to generate matrices and associated data for alphabets of > 2 symbols

	inputs:
	k_sym is the number of symbols (>2) (equal to the number of faces on the overall simplex)
	di_vals are the shift amounts corresponding to their respective symbol
	hk is the simplex height
	delta is the regular/random parameter, in range [0,1]
	p_and_s is 0 for pad only; 1 for pad and splice; 2 for splice only

	outputs:
	data_header, an Octave structure with fields,
		k_sym
		di_vals (vector of shift values [d_1,...,d_k])
		hk
		delta
		p_and_s
		nbrows (number of block rows)
		mxsz_vals (sizes of the matrices in block rows (vector of elements nbrows long))
	Si_mx_cell,
		a cell array of size [nbrows,k_sym]
		this stores, in sparse form, the k_sym Si binary matrices (entries in {0,1}) that constitute the fully regular sequence templates
		each cell row corresponds to a block row of the overall {S_1,...,S_k} matrix set; the row has k_sym matrices of [r,c] coordinates (of form [r,c;r,c;...]); column number j corresponds to symbol j
		block row i, if expanded from sparse form and into full binary matrices, will have matrices of dimension mxsz_vals(i)
	graph_array, an array of Octave graph structures, one structure for each block row, with fields,
		cyl_graph_matrix: a representation of the graph induced by the respective Si_mx_cell block row, in "cylindrical" wraparound format, with min_orbit steps
		binary_edge_matrix: alternate representation form for the graph, for use with fast isomorphism checking; stores the columnar binary string matrix, with entry/exit edges tallied in binary strings (as rows), of length 2*k_sym
		min_orbit: the number of steps used to make the graph cylindrical--ie think of modding it by min_orbit
		max_node_count: stores the maximum number of nodes over all steps in the post-min_orbit-modded graph form; this value is important for parsing the respective binary edge matrix; the following equality should always hold: (max node count)*(min_orbit)=(# of rows in respective binary edge matrix)
		node_count_vec: vector of length min_orbit, storing the node count at the respective time step
		node_numbering_vec: holds (columnar) vectors holding the original node numbers; it matches, row for row, the nodes referenced in the corresponding binary edge matrix; note, the binary edge graph will have the binary edge in/out strings ordered in (decreasing) lexicographical order, and the node_numbering_vec will match that ordering
		min_multiplicity_node: saves the node number of a node with binary edge code of minimal multiplicity for the graph
	Mi*P matrices (optional)
		as an option, the Mi*P matrices can be generated from the Si_mx_cell matrices, and saved in a .bin file
		they're stored by block rows, in row-major order--that is, the first block row, followed by the second block row, and so on
		the individual MiP matrix entries are written in row-major order as well (that is, (r1,c1),(r1,c2),...,(r2,c1),(r2,c2),...,...)

	nb: the indexing method to make matches requires the matrix of all possible coordinate k-tuples ("coord_mx") to be generated (here, by "sim_gen" recursive function) in decreasing lexicographical order

%}

clear global; % this seems necessary/important between runs of this function; some globals are used in sim_gen

global coord_mx; % for storing all possible k_sym-tuple altitude coordinates
global c_vals = zeros(1,k_sym);
global ix_cm = 1; % this will store one more than the number of coordinate k-tuples generated

BYTES_PER_INT16 = 2;

default_file_1 = ""; % optional; create default file name for saving S_i matrix set output in .m file
default_file_2 = ""; % optional; create default file name for saving M_i*P matrix set output in .bin file

if k_sym != length(di_vals)
	printf("problem with inputs\n");
	return
endif
if k_sym < 3
	printf("k_sym must be at least 3\n");
	return
endif



%***
% simplex coordinate generation
%***

hu = hk-1; % the height of the simplex in units of 1/(k-1) is entered as parameter hk; most of this routine uses "hu"--hu = h*(k-1)-1 = hk-1 (where h is the "original" simplex altitude)--eg for h=3 units, for a 2-simplex (3 symbols), hk=6, and hu=5; recall, this was for purposes of indexing from 0--ie to index along any simplex altitude by [0,1,2,...,hu], where 0 is the respective base of the simplex, and hu is at the base of the uppermost (apex) unit shear cell
tot_slots = hu+k_sym-1;

% notification of memory usage; the max memory used is by shft_coords cell array
% the number of weak compositions of n into k parts (sum to n of exactly k non-zero numbers) is choose(n+k-1,k-1) = choose(n+k-1,n); there should be (h(k-1)-1+k-1 choose k-1) = choose(hu+k_sym-1,k_sym-1) unit shear cells in all
num_initial_shift_cells = nchoosek(tot_slots,k_sym-1);
printf("%ld initial shift cell coordinates will require approx. %f GB\n",num_initial_shift_cells,k_sym*(num_initial_shift_cells*k_sym)*BYTES_PER_INT16/1024^3);

coord_mx = zeros(num_initial_shift_cells,k_sym,"int16"); % cast as int16's, to help reduce memory usage (the simplex height will be assumed to be less than 2^15-1)

% sim_gen will effectively list all sums of type a_1+...+a_{k_sym}=hu
sim_gen(hu,0,1,k_sym-1); % this gets: hu, initial offset (should be 0), starting recursion depth (should be 1), k-1 (1 less the recursion depth; this translates into the right level of recursion for k symbols)
% for a (k-1)-simplex of hk height units, this creates the ~symplectic k-tuples where each component can have values in [0,hu], representing the height of the unit shear cell from the corresponding face (k_sym total), and the "rule" is all the components in the k-tuple have to add to hu
% NB: sim_gen has to produce the coordinate k-tuples in decreasing lexicographic order (eg [5,0,0;4,1,0;4,0,1;3,2,0;3,1,1;3,0,2;...])--this is to exactly match the indexing produced by the (binomial coefficient) routine in the mx_shft post-shift coordinate matrix processing loop

% checks; if these trigger, something's wrong with the code
if num_initial_shift_cells!=ix_cm-1 % ix_cm came straight from the recursion routine sim_gen
	"wrong number of simplex coordinates"
	return;
endif
sum_chk = sum(coord_mx,2);
if sum(abs(sum_chk-hu))>1/num_initial_shift_cells
	"simplex coordinates not of constant sum"
	return;
endif



%***
% perform the symbol shift operations (stored in di_vals) on the simplex coordinate matrix, one for each symbol
%***

[rcm,ccm] = size(coord_mx); % coord_max has each row equal to a distinct k-tuple coordinate

shft_coords = cell(1,k_sym);
for j=1:k_sym
	shft_coords{j} = coord_mx; % the shift operators will each act separately on this initial coordinate set
endfor

psh_vec = ones(rcm,1,"int16"); % this is just [1,...,1]' (very long)
for j=1:k_sym
	shft_coords{j}(:,j) += di_vals(j)*(k_sym-1)*psh_vec; % ie this is how many units of 1/(k-1) we shift away from the face in question (the jth), which is applied to the jth column of coord_mx (copied into shft_coords)
	ix_rem = [1:k_sym];
	ix_rem(j) = [];
	for j1=1:length(ix_rem)
		shft_coords{j}(:,ix_rem(j1)) -= di_vals(j)*psh_vec; % all other faces have amount 1/(k-1) sheared off for each "step" (of size "1") in the shift amount (di_vals), ie shift value integer units in terms of coord_mx values
	endfor
endfor



%***
% initially process post-shifted coordinate matrices (k_sym in number, one for each symbol)
%***

% this saves shifted k-tuples that stay within the boundary simplex, and converts the post-shifted altitude coordinate k-tuples into cell index values, using an indexing scheme associated with binomial coefficients; this indexing was done for speed (saves having to search for matches)
% a row of the shft_coords matrices is rejected if it has any <0 or any >hu values
% the main result of this section is the mx_pairs cell array, each of k_sym entries a matrix of unit cell index pairs (of form [pre-shift index 1,post-shift index 1; pre-shift index 2, post-shift index 2;...])
mx_pairs=cell(1,k_sym);
for j=1:k_sym
	mx_shft = shft_coords{j}; % this is a matrix of row vectors, copied from coord_mx, after the (j-face) shift operation; each row contains the k_sym altitudes, each [0,hu]; the index of each row amounts to the original cell #, prior to the shift operation
	ct_mx_pr = 1;
	out_range = (mx_shft>hu) + (mx_shft<0);
	rw_rem = find(sum(out_range,2)>0.1); % column vector; which rows in shft_coords are no good
	ix_rr = [1:rcm]';
	ix_rr(rw_rem) = []; % this will leave a vector of indices in [1,#coord_mx rows] that match the post-row-removed shft_coords{j} matrix
	mx_shft(rw_rem,:) = [];
	[rmn,cmn] = size(mx_shft);

	% process the matrix mx_shft for this (jth) symbol; its rows are coordinate k-tuples that resulted from the shift of d_i on the simplex; need to determine the shear cell index (ie what row # of coord_mx) of the shifted k-tuples
	% this "thinks" in terms of binary equivalents to the coordinate k-tuples; eg [1,2,0,2] for k_sym=4 and hu=5 has binary-partition rep 01001100 (the "1"'s comparmentalize / partition the "0"'s)
	% to convert the binary rep to the index #, for the mth "1" placed, at position n, in from left (starting at 1), we add nchoosek(tot_slots-n,(k_sym-1)-m+1) (sum over all m in [1,k_sym-1])
	% nb: it bears repeating this process requires sym_gen enumerate the rows in coord_mx in a specific way, with values down a column always decreasing as row # increases--what amounts to strict decreasing lexicographical order (eg [9,0,0;8,1,0;8,0,1;7,2,0;7,1,1;7,0,2;...])--if this weren't the case, could just run this binomial indexing process on both coord_mx, and on the individual (k_sym in all) mx_shft matrices

	tmp_mx = mx_shft(:,1:end-1); % we only want columns 1:k_sym-1
	posn_n = cumsum(tmp_mx+1,2); % now the rows represent the indices of the "1"'s in the binary-partition rep
	num_rows = length(tmp_mx(:,1));

	mth_1 = ones(num_rows,1)*[1:k_sym-1]; % at column m, we've placed the mth "1"
	tot_ones = k_sym-1; % how many 1's there are to distribute (in binary partition format 010001010, eg; see HP draft 1 supp looseleaf p. 81)

	% create two matrices of size (num_rows,tot_ones), to hold the values for the n-choose-k operation (can think of this as setting up an Octave "dot" form of n-choose-k, as with [matrix 1]./[matrix 2])
	nc_upper = tot_slots-posn_n; % upper value of choose((hu+k_sym-1)-n, (k_sym-1)-m+1)
	nc_lower = tot_ones-mth_1+1; % denominator of choose((hu+k_sym-1)-n, (k_sym-1)-m+1); this should always just be a matrix of num_rows identical rows, [k_sym-1,k_sym-2,...,1]

	null_ix = ((nc_upper-nc_lower)<0); % index all nc_upper,nc_lower value pairs that a-choose-b can't be performed on (ie those where b is greater than a)

	nc_upper(null_ix) = nc_lower(null_ix); % this makes the next steps easier, setting all to-be-nulled upper/lower pairs to (c,c), where "c" is the column number; these values (nchoosek(c,c)=1) still have to be eliminated from the row sums

	nck = [];
	for tt=1:tot_ones % go through the columns of nc_upper and nc_lower, performing the n-choose-k operation

		col_upper = nc_upper(:,tt); % pull tt'th column from col_upper (the "a"'s of a-choose-b)
		tmp_max_sub = tot_ones-tt+1; % for creating the vector a-[0:tmp_max_sub-1] out of each "a" of a-choose-b
		tmp_1 = ones(num_rows,1)*[0:tmp_max_sub-1]; % each row has the vector [0:tmp_max_sub-1]
		tmp_2 = col_upper*ones(1,tmp_max_sub);
		tmp_upper = tmp_2-tmp_1; % each row has the vector a-[0:tmp_max_sub-1], where a is the entry from col_upper
		tmp_lower = nc_lower(:,tt:end); % this is the b-[0:tmp_max_sub-1] vector, for a-choose-b
		nck(:,tt) = round(prod(tmp_upper./tmp_lower,2)); % this takes the products of vectors ((a-[0:tmp_max_sub-1])./(b-[0:tmp_max_sub-1])); for sufficiently large values, this can run into precision errors

	endfor

	nck(null_ix) = 0; % this zeros the null entries
	nck_ix_out = 1+sum(nck,2); % this amounts to assigning the rows of coord_mx that were not shifted to null under the d_i shift (those rows indexed by column vector ix_rr) their proper shear cell indices
	mx_pairs{j} = [nck_ix_out,ix_rr];

endfor



%***
% check all shifts have mapped at least one unit cell into another
%***

tmp_ptp = 0;
for jj=1:k_sym
	if isempty(mx_pairs{jj})
		printf("initial shift amounts have mapped one or more entire simplices to null; try a larger simplex height\n");
		return
	endif
	tmp_ptp += size(mx_pairs{jj},1);
endfor

printf("%ld total shift pairs to initially process\n",tmp_ptp);



%***
% pare all columns of zeros, and all rows of zeros
%***

post_pare = MultAlpha_Pare(num_initial_shift_cells,mx_pairs); % post_pare will be a cell array row vector of size k_sym+1; the first cell has a single number, the matrix dimension after the paring step (it's set to -1 if there are no shift cells remaining after the paring step); cells 2 thru k_sym+1 hold the respective symbols' matrices of coordinate pairs (format, [r1,c1;r2,c2;...])

if (post_pare{1}==-1)
    printf("nothing left after paring step; try a larger simplex height\n");
    return
endif

mx_dim = post_pare{1};
tmp_ptp = 0;
routes_in = cell(1,k_sym);
for j=1:k_sym
	routes_in{j} = post_pare{j+1};
	tmp_ptp += size(routes_in{j},1);
endfor

printf("%ld total shift pairs to process after paring step\n",tmp_ptp);



%***
% find all closed routes in the k_sym binary matrices of routes_in
%***

[Si_mx_cell,graph_array] = MultAlpha_Routes(mx_dim,routes_in,di_vals); % finds any/all closed cycles in the k_sym binary matrices recorded in sparse form in routes_in cell array; produces the Si_mx_cell array, portioned into closed-route block rows, and the graph_array array of graph structures



%***
% display info from the Si_mx_cell's block rows (optional); parse Si_mx_cell into (sparse-format) matrices
%***

[om_r,om_c] = size(Si_mx_cell); % each row corresponds to a given row in the main shift matrices' (S_i) block diagonals, and has k_sym columns

if om_r==0 % notification the simplex height (hk) has been chosen too small for this shift set (there are no closed routes that "fit" in the simplex)
	"hk value too small for given d_i"
	return
endif

printf("%i block rows\n",om_r);

for j1=1:om_r
	mpr = 0;
	for j2=1:om_c
		tmp_max = max(max(Si_mx_cell{j1,j2})); % finds max of entries of sparse matrix form [r1,c1;r2,c2;...] for all symbols in this block row; this is the matrix dimension for the respective (j1) block row
		if tmp_max>mpr
			mpr = tmp_max;
		endif
	endfor
	mpr_vec(j1) = mpr;
endfor

mpr_unq = unique(mpr_vec);

num_entries = 0;
for j1=1:length(mpr_unq)
	tmp_ct = sum(mpr_vec==mpr_unq(j1));
	printf("%i of size %i x %i\n",tmp_ct,mpr_unq(j1),mpr_unq(j1));
	num_entries += tmp_ct*(k_sym*(mpr_unq(j1))^2); % totals number of entries (elements) in all the matrices of all the block rows
endfor



% ***
% save file; optional user input--can be commented out, hard coding fname_1 a/o fname_2 in this file and saving automatically
% ***

fflush (stdout);
str_yn = input("Save the results (either the S_i, in a compact .m file, or the full M_i*P forms, in a .bin file) (y/n)? ","s");
if str_yn!="Y" && str_yn!="y"
	return;
endif



%***
% option to save the S_i in sparse form to an .m file
%***

str_yn = input("Save the S_i and graph_array to an Octave .m file (y/n)? ","s");
if str_yn=="Y" || str_yn=="y"
	if !strcmp(default_file_1,"")
		printf("Enter filename (return for default file, %s): ",default_file_1);
		fname_1 = input("","s");
		if strcmp(fname_1,"")
			fname_1 = default_file_1;
		endif
	else
		fname_1 = input("Enter filename: ","s");
	endif
	data_header.k_sym = k_sym;
	data_header.di_vals = di_vals;
	data_header.hk = hk;
	data_header.delta = delta;
	data_header.p_and_s = p_and_s;
	data_header.br_dims = mpr_vec;
	save(fname_1,'data_header','Si_mx_cell','graph_array');
endif



%***
% option to save the M_i*P in full form to a .bin file
%***

% because this fully unpacks the block rows out of sparse form, this can consume a lot of memory
bytes_per_double = 8;
MB_req = num_entries*bytes_per_double/(1024^2);
printf("Save the M_i*P products to a bin file (file will require approximately %.2f MB) (y/n)? ", MB_req);
str_yn = input("","s");
if str_yn!="Y" && str_yn!="y"
	return;
endif
if !strcmp(default_file_2,"")
	printf("Enter filename (return for default file, %s): ",default_file_2);
	fname_2 = input("","s");
	if strcmp(fname_2,"")
		fname_2 = default_file_2;
	endif
else
	fname_2 = input("Enter filename: ","s");
endif



% ***
% unpack S_i sparse matrix cell array into M_i*P (block diagonal) products
% ***

fid = fopen(fname_2,'w');

	% write the header info
fwrite(fid,k_sym,'double'); % <k_sym> number of symbols
fwrite(fid,di_vals,'double'); % <di_vals> (vector) writes shift amounts
fwrite(fid,hk,'double'); % <hk> writes height of simplex
fwrite(fid,delta,'double'); % <delta>
fwrite(fid,p_and_s,'double'); % =0 for pad only; =1 for pad and splice; =2 for splice only
fwrite(fid,om_r,'double'); % <nbrows> number of block rows
fwrite(fid,mpr_vec,'double'); % <mxsz_vals> (vector) sizes of the matrices in the block rows

	% write the Mi*P matrix products
p_vals = (1./di_vals)/(sum(1./di_vals)); % probabilities for symbol #s 1,2,...,k
si_br_cells = cell(1,k_sym);
for j=1:om_r % go through each block row of the S_i
	mx_sz = mpr_vec(j); % this was the sizes of the matrices in the block row
	tmp_sum_mx = zeros(mx_sz);
	for k=1:om_c % go through each of the k_sym matrices in the block row, reconstituting the S_i block diagonal matrices into full (non-sparse) form
		tot_en = size(Si_mx_cell{j,k},1);
		si_br_cells{k} = full(sparse(Si_mx_cell{j,k}(:,1),Si_mx_cell{j,k}(:,2),ones(1,tot_en),mx_sz,mx_sz));
		tmp_sum_mx += delta*p_vals(k)*si_br_cells{k}; % for creating the padding matrix (if p_and_s != 2)
	endfor
	if (p_and_s==0 || p_and_s==1)
		pad_mx = inv(eye(mx_sz)-tmp_sum_mx);
	endif
	for k=1:om_c
		if p_and_s==0 % pad only
			MiP = si_br_cells{k}*pad_mx;
		elseif p_and_s==1 % pad and splice
			MiP = (delta*p_vals(k)*eye(mx_sz) + (1-delta)*si_br_cells{k}) * pad_mx;
		else % splice only
			MiP = delta*p_vals(k)*eye(mx_sz) + (1-delta)*si_br_cells{k};
		endif
		mx_wrt = conj(MiP'); % this to render matrix in row-major form in the bin file
		fwrite(fid,mx_wrt,'double');
	endfor
endfor

fclose(fid);

endfunction







function sim_gen(hu,x_1,rd,rd_max)

	% this function is called recursively, to generate all k-tuples {a_1,...,a_k} such that a_1+...+a_k=hu, which provides the altitude coordinates of every unit shear cell in the simplex of height hu+1 units

	global coord_mx;
	global c_vals;
	global ix_cm;

	for x_2=(hu-x_1):-1:0
		if rd<rd_max
			c_vals(rd) = x_2;
			sim_gen(hu,x_1+x_2,rd+1,rd_max);
		else
			coord_mx(ix_cm,:) = [c_vals(1:rd_max-1),x_2,hu-(x_1+x_2)];
			ix_cm++;
		endif
	endfor

endfunction






