function BinAlpha_Le_Var_MC(prod_len,num_sam,MiP_bin_file)

%{

	a Monte Carlo matrix product routine, binary alphabets, for estimating the Lyapunov exponent and variance

	inputs:
	prod_len is the number of matrices in the product chains
	num_sam is the number of product chains
	MiP_Bin_file is a string, the name of the .bin file holding the {M_1*P,M_2*P} matrix product pair

	outputs:
	produces LE and var estimates, with option to save the [LE est, var est, prod_len, num_sam] vector in an .m file

	this uses the LV_MC_BA C++ routine, called via a system call
	nb: as of latest, this requires the called C++ routine to produce as its first two lines of numeric output (strings / non-numeric lines notwithstanding), the LE est and the var est, respectively; eg the C++ routine output could be "LE est:\n1.235\nvar est:\n3.212", but not "LE est: 1.235\nvar est: 3.212"
	nb: this version requires the user hard-coding the name/location of the LV_MC_BA.cpp executable

%}


%***
% for user-supplied C++ Monte Carlo file name and location
%***
LV_MC_BA_file = ""; % non-optional, for user to indicate the C++ MC routine file (generated from LV_MC_BA.cpp)
if isnull(LV_MC_BA_FILE)
	printf("assign filename string to LV_MC_BA_FILE in code for LV_MC_BA.cpp executable\n");
	return
endif

default_file = ""; % optional, for default save output file


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
	br_dims(jj) = fread(fid,1,'double'); % <mxsz_vals> (vector) sizes of the matrices in the block rows; this should always just be one element, equal to hk
endfor
fclose(fid);

 	% (optional) display .bin file header info:
%{
printf("%i symbols; hk=%i; delta=%f; p_and_s=%i; %i block rows; ",...
	k_sym,hk,delta,p_and_s,nbrows);
printf("di_vals={");
for jj=1:k_sym-1
	printf("%i,",di_vals(jj));
endfor
printf("%i}\n",di_vals(k_sym));
%}



%***
% call the C++ Monte Carlo routine, via system call
%***
str_in = [LV_MC_BA_file " " num2str(prod_len) " " num2str(num_sam) " " MiP_bin_file];
[status, output] = system(str_in);


%***
% parse the output string, produced from the system call to the C++ routine
%***

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
if (ct<2)
	printf("error in C++ Monte Carlo routine output; must display LE and var estimates on their own lines\n");
	return
endif

fam_LE_var_mx = [LE_est,var_est,prod_len,num_sam]; % the variable naming (fam_LE_var_mx) is for consistency with the same LE / var routine for alphabets with > 2 symbols


%***
% display results and offer option to save
%***

printf("LE: %f\n var: %f\n",LE_est,var_est);

fflush(stdout);

printf("\n");

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







