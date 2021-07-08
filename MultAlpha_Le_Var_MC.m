function MultAlpha_Le_Var_MC(arg1, arg2, arg3, arg4, arg5)

%{

	Monte Carlo estimator for the Lyapunov exponent and variance, for alphabets of more than 2 symbols
	the routine can produce estimates for one or more isomorphism families, selecting an arbitrary block row within the given family

	inputs:
	arg1--product length; number of matrices in Monte Carlo products
	arg2--number of samples; number of random matrix products in the Monte Carlo sampling
	arg3--Mi*P matrices .bin filename
	arg4--isomorphism families .m filename
	arg5--(optional) selection vector of isomorphism families; if this argument is omitted, the routine will produce estimates for all block row isomorphism families in br_iso_families
		the selection vector can contain a list of values corresponding to which isomorphism families the Lyapunov exponent and variance are to be estimated for
		values in the selection vector correspond to (outer) indices in the br_iso_families cell array, which enumerates isomorphic block row families, in order of increasing block row matrix dimension

	output:
	fam_LE_var_mx--a matrix, with rows corresponding to isomorphism families in the order they appear in the br_iso_families cell array; each row holds the Lyapunov exponent estimate, the variance estimate, the product length used, and the number of samples used; if sel_vec is inputted, the families not used will have a 0 on all 4 columns

%}

prod_len = arg1;
num_sam = arg2;
MiP_bin_file = arg3;
isom_cell_file = arg4;

if (nargin==5)
	sel_vec = arg5;
else
	sel_vec = [];
endif

default_file = ""; % optional, for default save output file



%***
% for user-supplied C++ Monte Carlo file name and location
%***

LV_MC_MA_file = ""; % non-optional, for user to indicate the C++ MC routine file (generated from LV_MC_BA.cpp)
if isnull(LV_MC_MA_FILE)
	printf("assign filename string to LV_MC_MA_FILE in code for LV_MC_MA.cpp executable\n");
	return
endif



%***
% open isomorphism family .m file, and Mi*P matrix .bin files
%***

try
	load(isom_cell_file,"data_header","br_iso_families");
catch
	printf("problem with input filename, isomorphism family .m file\n");
	return
end_try_catch

printf("\n");



%***
% open the Mi*P binary file, and retrieve the header info
%***

fid = fopen(MiP_bin_file);
if (fid==-1)
	printf("problem with Mi*P .bin file\n");
	return
endif
k_sym = fread(fid,1,'double'); % <k_sym> number of symbols
for jj=1:k_sym
	di_vals(jj) = fread(fid,1,'double'); % <di_vals> (vector) shift amounts
endfor
hk = fread(fid,1,'double'); % <hk> height of simplex
delta = fread(fid,1,'double'); % <delta>
p_and_s = fread(fid,1,'double'); % =0 for pad only; =1 for pad and splice; =2 for splice only
nbrows = fread(fid,1,'double'); % <nbrows> number of block rows
for jj=1:nbrows
	br_dims(jj) = fread(fid,1,'double'); % <mxsz_vals> (vector) sizes of the matrices in the block rows
endfor
fclose(fid);

	% optional check; check the .bin Mi*P header info against the data_header from the isom_cell_file
if !(data_header.k_sym==k_sym)||!isequal(data_header.di_vals,di_vals)||!(data_header.hk==hk)||!(data_header.delta==delta)||!(data_header.p_and_s==p_and_s)||!isequal(data_header.br_dims,br_dims)
	printf("file data headers don't agree between .m and .bin files\n")
	return
else
	printf("%i symbols; hk=%i; delta=%f; p_and_s=%i; %i block rows; ",...
		data_header.k_sym,data_header.hk,data_header.delta,data_header.p_and_s,length(data_header.br_dims));
	printf("di_vals={");
	for jj=1:k_sym-1
		printf("%i,",data_header.di_vals(jj));
	endfor
	printf("%i}\n",data_header.di_vals(k_sym));
endif



%***
% read in info from br_iso_families, and determine which block row #s will be used in the main Monte Carlo loop
%***

sz_icf = size(br_iso_families);
if (min(sz_icf)!=1)
	printf("br_iso_families file should hold a vector cell array\n");
	return
endif
num_mx_dims = max(sz_icf); % how many different-sized matrices are represented in br_iso_families

% load a representative block row index into "iso_fam_rep" from each family in the br_iso_families inner cell arrays
ct = 0;
for jj=1:num_mx_dims
	cell_tmp = br_iso_families{jj};
	sz_tmp = size(cell_tmp); % individual entries in the br_iso_families cell vector are themselves (vector) cell arrays
	if (min(sz_tmp)!=1)
		printf("br_iso_families cell array vector members should hold a(nother) vector cell array\n");
		return
	endif
	for kk=1:max(sz_tmp)
		ct++;
		iso_fam_rep(ct) = cell_tmp{kk}(1); % choose an arbitrary element (block row index) from this family (may as well be the first element in the cell_tmp{kk} vector of values)
	endfor
endfor

if !isempty(sel_vec)
	if (max(sel_vec)>length(iso_fam_rep) || min(sel_vec)<1)
		printf("problem with sel_vec input\n");
		return
	endif
	br_ix = iso_fam_rep(sel_vec);
else
	br_ix = iso_fam_rep;
endif



%***
% main loop over relevant isomorphic family block rows
%***

fam_LE_var_mx = zeros(length(iso_fam_rep),4);

for jj=1:length(br_ix)

	% call the Monte Carlo C++ routine, via system call
	str_in = [LV_MC_MA_file " " num2str(prod_len) " " num2str(num_sam) " " MiP_bin_file " " num2str(br_ix(jj))];
	[status, output] = system(str_in);

	% parse the output string, produced from the system call to the C++ routine
	% this looks for the first two lines (delimited by '\n') that have numbers in them--taking the 1st such to be the estimate of the Lyapunov exponent, and the second to be the variance estimate
	out_parse = strsplit(output,'\n'); % output will be one long string, with '\n' inserted at line breaks; we want the two lines with numbers--those with the LE and var estimates
	ct = 0;
	for kk=1:length(out_parse)
		tmp_val = str2num(out_parse{kk});
		if !isempty(tmp_val)
			ct++;
			if (ct==1)
				LE_est = tmp_val;
			else
				var_est = tmp_val;
			endif
		endif
		if (ct>=2)
			break
		endif
	endfor
	% nb: as of latest, this requires the C++ routine to produce as its first two lines of numeric output (strings / non-numeric lines notwithstanding), the LE est and the var est, respectively; eg the C++ routine output could be "LE est:\n1.235\nvar est:\n3.212", but not "LE est: 1.235\nvar est: 3.212"
	if (ct<2)
		printf("error in C++ Monte Carlo routine output; must display LE and var estimates on their own lines\n");
		return
	endif

	if !isempty(sel_vec)
		printf("family #%i, LE: %f, var: %f\n",sel_vec(jj),LE_est,var_est);
		fam_LE_var_mx(sel_vec(jj),1:4) = [LE_est,var_est,prod_len,num_sam];
	else
		printf("family #%i, LE: %f, var: %f\n",jj,LE_est,var_est);
		fam_LE_var_mx(jj,1:4) = [LE_est,var_est,prod_len,num_sam];
	endif

	fflush(stdout);

endfor

printf("\n");


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

save(fname,'data_header','fam_LE_var_mx');



endfunction







