function BinAlpha_Mx_Gen_1(mx_sz,d_1,d_2,delta,p_and_s)

%{

	generates Mi*P matrices for binary alphabets ({0,1})

	inputs:
		mx_sz--proportion spread (analagous to simplex height in alphabets > 2 symbols)
		d_1, d_2--shift amounts (with d_1<=d_2, p_1=P(0)=d_2/(d_1+d_2)>=1/2, and d_1+d_2<=mx_sz)
		delta--regularity/randomness parameter, (in [0,1])
		p_and_s--pad and splice flag (p_and_s: 0 is pad only; 1 is pad-and-splice; 2 is splice only)

	outputs:
	this has the option to save the resulting {M_1*P,M_2*P} matrix product pair in a .bin file, 
	.bin file has standard header:
	<k_sym>, <di_vals>, <hk>, <delta>, <p_and_s>, <nbrows>, <mxsz_vals>
	the header is followed by the matrix entries, in row-major order (that is, (r1,c1),(r1,c2),...; (r2,c1),(r2,c2),...;...)
	nb: for binary alphabets, k_sym=2, di_vals=[d_1,d_2], hk=mx_sz, nbrows=1, and mxsz_vals is a single element vector, [mx_sz]
	nb: all values (header and matrix entries) are recorded in double precision (8 byte standard, 8 bits to a byte, on Octave)
		in general, for matrix data file purposes, this expects
		(1) bytes to be 8 bits
		(2) float type "double" to be consistent between here (the Octave routine) and the data file-receiving C++ routine(s) (otherwise the C++ routine will not be able to read the passfile properly)
		(there are some checks for this in the code of the relevant C++ files)

%}



default_file = ""; % optional; can set a default file here for saving


if (d_1>d_2) % reverses these, so that d_1<=d_2, just in case
	tmp = d_1;
	d_1 = d_2;
	d_2 = tmp;
endif
p_vals(1) = d_2/(d_1+d_2);
p_vals(2) = 1-p_vals(1);
k_sym = 2;
if (mx_sz < d_1+d_2)
	"ERROR in mx_sz, d_1, d_2 values"
endif


% ****
% create M_i P matrix products
% ****
S_0 = diag(ones(1,mx_sz-d_1),d_1);
S_1 = diag(ones(1,mx_sz-d_2),-d_2);
p = p_vals(1);
q = p_vals(2);
if p_and_s == 0	% pad only (P S_i P S_i ... P S_i P, w/ P = I + delta(pS_0+qS_1) + delta^2 (pS_0+qS_1)^ + ...)
	M_0 = S_0; % ie M_0 = S_0 under pad only
	M_1 = S_1; % ie M_1 = S_1 under pad only
	pad_0 = delta*p*S_0;
	pad_1 = delta*q*S_1;
	P = inv(eye(mx_sz)-(pad_0+pad_1));
elseif p_and_s == 1 % pad and splice (P M_i P M_i ... P M_i P, w/ M_i = p_i*delta*I + (1-delta)S_i, and P = I + delta(pS_0+qS_1) + delta^2 (pS_0+qS_1)^ + ...)
	M_0 = p*delta*eye(mx_sz)+(1-delta)*S_0; % eye(mx_sz) component allows splicing
	M_1 = q*delta*eye(mx_sz)+(1-delta)*S_1; % eye(mx_sz) component allows splicing
	pad_0 = delta*p*S_0;
	pad_1 = delta*q*S_1;
	P = inv(eye(mx_sz)-(pad_0+pad_1));
elseif p_and_s == 2 % splice only (P M_i P M_i ... P M_i P, w/ M_i = p_i*delta*I + (1-delta)S_i, and P=I)
	M_0 = p*delta*eye(mx_sz)+(1-delta)*S_0; % eye(mx_sz) component allows splicing
	M_1 = q*delta*eye(mx_sz)+(1-delta)*S_1; % eye(mx_sz) component allows splicing
	P = eye(mx_sz); % no padding for splice only
endif
	% these are the matrix product pairs passed to the C++ routine via the passfile (their transposes are taken before writing, to retain read-by-row (left-right -> up-down) order; Octave uses column-major order to read/write, while in the C++ code it's more natural to use row-major order for the matrices)
MP_i_vals = cell(1,2);
MP_i_vals{1} = M_0*P;
MP_i_vals{2} = M_1*P;



% ****
% option to save the Mi*P matrix pair
% ****
printf("Save the M_i*P products to a bin file (y/n)? ");
str_yn = input("","s");
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


% ****
% write to .bin file
% ****
fid = fopen(fname,'w');

	% write header
fwrite(fid,2.0,'double'); % <k_sym> writes number of symbols
fwrite(fid,[d_1,d_2],'double'); % <di_vals> writes shift amounts
fwrite(fid,mx_sz,'double'); % <hk> writes height of "simplex"; in binary sequences this is always the size of the involved matrices
fwrite(fid,delta,'double'); % <delta>
fwrite(fid,p_and_s,'double'); % <p_and_s>
fwrite(fid,1,'double'); % <nbrows> number of block rows; for binary sequences this is always just 1
fwrite(fid,mx_sz,'double'); % <sz_1,...,sz_t> sizes of the block rows; for binary sequences this again is just the matrix size

	% write matrices
for i=1:2
	wrt_prd_i = conj(MP_i_vals{i}'); % conjugate-transpose for recording in row-major order
	for j=1:numel(wrt_prd_i)
		count = fwrite(fid,wrt_prd_i(j),'double'); % octave vertically "vectorizes" an array by [1,3;2,4]--ie runs down the 1st column, then 2nd, etc (top to bottom, left to right)--so to get right-left->up-down have taken the transpose first
	endfor
endfor

fclose(fid);




endfunction
