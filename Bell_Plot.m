function Bell_Plot(prod_len,num_sam,lambda,var,MiP_bin_file)

%{

	compares a random sample of matrix products to the theoretical bell curve derived from Lyapunov exponent and variance estimates (from central limit theorem)
	works for any alphabet
	after running through all matrix product samples, a normalized frequency plot of the (logs of) norms is made along the log(||.||) axis; the theoretical normal curve (mean=prod_len*lambda; variance=prod_len*var) is also plotted

	inputs:
	prod_len--the length of the matrix products under random sampling
	num_sam--how many matrix products to compute
	lambda--estimated Lyapunov exponent for the Mi*P matrix set
	var--estimated variance for the Mi*P matrix set
	MiP_bin_file--name of .bin file holding the Mi*P matrix set

	outputs:
	adj_factor--this is prod_len*(sample mean - lambda); this allows estimating an adjustment for the better computation of z-scores, especially for shorter matrix product lengths, near the rule-of-thumb lower limit of 10 times the largest block row's matrix dimension; given log(||.||), the log of the entrywise norm of the matrix product being measured, the z-score estimate becomes:
		(log(||.||)-(prod_len*lambda+adj_factor)) / sqrt(prod_len*var)
	plot showing both the frequence histogram and the normal distribution

%}


%***
% open the MiP_bin_file binary file, and read in the header info, and the MiP matrices
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

p_vals = (1./di_vals)/(sum(1./di_vals)); % probabilities for symbol #s 1,2,...,k
p_cdf = cumsum(p_vals);
p_cdf(end) = 1.0; % safeguard against precision errors



% ****
% Monte Carlo loop through num_sam
% ****

for kk=1:num_sam

	% optional, but helps reduce storage
	vec_cell = cell(nbrows,1);
	for jj=1:nbrows
		vec_cell{jj} = ones(mxsz_vals(jj),1);
	endfor

	log_vec_norm = zeros(nbrows,1);

	for mm=1:prod_len

		sym_now = min(find(rand()<p_cdf)); % generalized Bernoulli random selection of the next symbol

		for jj=1:nbrows % go through each block row
			vec_cell{jj} = tot_MiP_mx_cell{jj,sym_now}*vec_cell{jj};
			tmp_nrm = norm(vec_cell{jj});
			vec_cell{jj} = vec_cell{jj}/tmp_nrm;
			log_vec_norm(jj) += log(tmp_nrm);
		endfor % end loop over block rows

	endfor % end loop over matrix product

	% find the log of the entrywise norm for the whole matrix (all block rows) from individual block row log(||.||) values
	for jj=1:nbrows
		tmp_log_norm(jj) = log_vec_norm(jj)+log(sum(vec_cell{jj})); % vector of log(||.||)'s, for each block row's product
	endfor
	tmp_offset = max(tmp_log_norm);
	tot_log_norm(kk) = tmp_offset + log(sum(exp(tmp_log_norm-tmp_offset)));

endfor % end loop over num_sam



%***
% compute and display adjustment factor
%***

sam_mean = mean(tot_log_norm)/prod_len; %ie this is the mean of all our sample values; this gives the approximation to the Lyapunov exponent
adj_factor = prod_len*sam_mean - prod_len*lambda;

printf("theoretical mean: %f; sample mean: %f; adjustment factor est: %f\n", prod_len*lambda, prod_len*sam_mean, adj_factor);



%***
% create and display histogram of random matrix product samples
%***

num_hist_bins = 40;
[num_el, bins] = hist(tot_log_norm,num_hist_bins); % num_el will contain the frequencies (of the norms) within respective ranges (bins)
tot_el = sum(num_el); % should(?) equal num_sam

% to normalize, divide num_el by respective bin width
tmp_a = (bins(2)-bins(1))/2;
tmp_b = (bins(end)-bins(end-1))/2;
tmp_c = bins(2:end)-bins(1:end-1);
bin_edges = bins(1:end-1)+tmp_c/2;
bin_edges = [bins(1)-tmp_a bin_edges bins(end)+tmp_b];
num_norm = num_el ./ (bin_edges(2:end)-bin_edges(1:end-1));

figure;
plot(bins,num_norm/tot_el,'k'); % normalized histogram



%***
% compute and show normal distribution corresponding to lambda and var
%***

mesh = sqrt(prod_len*var)/10;
bins = [-4*sqrt(prod_len*var)+prod_len*lambda:mesh:4*sqrt(prod_len*var)+prod_len*lambda];
hist_nor = exp(-(bins-prod_len*lambda).^2/(2*var*prod_len))/sqrt(prod_len*2*pi*var);
hold on;
plot(bins,hist_nor,'r');


endfunction















