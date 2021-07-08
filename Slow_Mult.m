function Slow_Mult(sym_string,MiP_bin_file)
%{

	uses Octave to calculate matrix product of the input symbol string; matrices are read from a standard-format .bin file
	for a binary alphabet, the symbols must be in {0,1}
	for alphabets of more than two symbols, symbols must be in {1,...,k_sym}

	inputs:
	sym_string--a string of symbols, to be rendered as a matrix product
	MiP_bin_file--name of .bin file holding the set of Mi*P matrix product pairs

	output:
	tot_log_norm--log of the entrywise norm of the matrix product corresponding to the input sym_string

%}


%***
% open the Mi*P .bin file, and retrieve the header info and Mi*P set
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
	mxsz_vals(jj) = fread(fid,1,'double'); % <mxsz_vals> (vector) sizes of the matrices in the block rows
endfor
for jj=1:nbrows
	for kk=1:k_sym
		tot_MiP_mx_cell{jj,kk} = (fread(fid,[mxsz_vals(jj),mxsz_vals(jj)],'double'))'; % to unpack the matrices from the bin file, recall Octave uses column-major form, while the matrices have been stored in row-major form; so transpose the matrix read w/ fread before saving in tot_MiP cell
	endfor
endfor
fclose(fid);



%***
% initializations
%***

% optional, but helps reduce storage
vec_cell = cell(nbrows,1);
for jj=1:nbrows
	vec_cell{jj} = ones(mxsz_vals(jj),1);
endfor

log_vec_norm = zeros(nbrows,1);
% eg input symbols "1233" reads M_1*P*M_2*P*M_3*P*M_3*P; since this follows the literature in building matrix products out "from the right" (each additional matrix drawn multiplies from the left; ie S_m*...*S_1), to compute the product we need to draw the matrices off the string in "reverse" order (from right to left), multiplying the accumulation vector, in turn, by each
if (k_sym==2) % for binary alphabets
	sym_vec = fliplr(sym_string-'0'+1);
else
	sym_vec = fliplr(sym_string-'1'+1);
endif



%***
% matrix multiplication loop
%***

for mm=1:length(sym_vec) % loop over symbols

	sym_now = sym_vec(mm);
	if (sym_now<1 || sym_now>k_sym)
		printf("problem with input symbol string; symbols out of range\n");
		return
	endif

	for jj=1:nbrows % go through each block row
		vec_cell{jj} = tot_MiP_mx_cell{jj,sym_now}*vec_cell{jj};
		tmp_nrm = norm(vec_cell{jj});
		vec_cell{jj} = vec_cell{jj}/tmp_nrm;
		log_vec_norm(jj) += log(tmp_nrm);

	endfor % end loop over block rows

endfor % end loop over sym_string



%***
% output log of entrywise norm of product, and exit
%***

% find the log of the entrywise norm for the whole matrix (all block rows) from individual block row log(||.||) values
for jj=1:nbrows
	tmp_log_norm(jj) = log_vec_norm(jj)+log(sum(vec_cell{jj})); % vector of log(||.||)'s, for each block row's product
endfor
tmp_offset = max(tmp_log_norm);
tot_log_norm = tmp_offset + log(sum(exp(tmp_log_norm-tmp_offset)));

printf("log||.||=%f\n",tot_log_norm);



endfunction
