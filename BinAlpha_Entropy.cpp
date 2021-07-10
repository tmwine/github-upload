/*

	entropy estimates for binary alphabets, using Vanneste's resampled Monte Carlo algorithm

	inputs (command line):
		prod_len--number of matrices (symbols) in product
		num_sam--number of sample products
		<Mi*P data filename>--the appropriate .bin file containing the relevant {M0*P,M1*P,P} matrix set, with standard data header

	output:
		entropy estimate

	this reads in a standard format M_i*P matrix data file:
	<k_sym>, <di_vals>, <hk>, <delta>, <p_and_s>, <nbrows>, <mxsz_vals>, <MiP>
	k_sym = number of symbols
	di_vals = shift values (d_1,...,d_k)
	hk = simplex height (=matrix size in binary alphabet case)
	delta = pad/splice penalty parameter, for degree of regularity/randomness 
	p_and_s = 0 for pad only, 1 for pad and splice, 2 for splice only
	nbrows = number of block rows (=1 for binary alphabets)
	mxsz_vals = sizes of the matrices in respective block rows (array of elements nbrows long)
	MiP = product matrix data; it's stored in order MiP_1,1,...,MiP_1,k; MiP_2,1,...,MiP_2,k; ...; ie the first block row, followed by the second block row, and so on
	(P = padding matrix data)
	nb: the individual MiP matrix entries are expected to be written in row-major order (that is, (r1,c1),(r1,c2),...,(r2,c1),(r2,c2),...,...)

	nb: this should compile on C++11 or later; eg $g++ -std=c++11 my_code.cpp -o my_code

*/


using namespace std; // this is generally bad practice

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string.h>
#include <sys/time.h>
#include <iomanip> // for setprecision (at least)
#include <climits> // for checking bits in a byte
#include <fstream> // for string file access

#define PAD_ONLY 0
#define PAD_AND_SPLICE 1
#define SPLICE_ONLY 2
#define LINES_PER_BLOCK 3 // this determines how many file lines to read before performing the matrix multiplication (ie processing in batches)


void cout_vec(double *vec, int s); // for debugging


int main(int argc,char *argv[])
{

	if (argc != 4) {
		cout << "Problem with input arguments." << endl;
		exit(1);
	}

	long prod_len = atol(argv[1]);
	long num_sam = atol(argv[2]);
	char *MiPfile = argv[3];

	// check C++ here is using 8 bit bytes, and 8 bytes to a "double" (for Octave .bin file compatability)
	if (CHAR_BIT!=8) {
		cout << "Number of bits in a byte here is " << CHAR_BIT << ", not 8. This will likely produce problems with reading the Octave-generated M_i*P .bin data file, which renders doubles with 64 bits. Exiting." << endl;
		exit(1);
	} else if (sizeof(double)!=8) {
		cout << "C++ defines type double as " << sizeof(double) << " bytes. This will likely produce problems reading the Octave-generated M_i*P .bin data file, which renders doubles with 64 bits. Exiting." << endl;
	}



	//***
	// read in .bin file header
	//***

	// get prelim parameters for matrix set (for binary alphabets, k_sym=2, nbrows=1)
	// the matrix data file has (in order) {matrices' size, number of symbols, p_1,...,p_k, M_1,...,M_k}

	FILE *pfile = NULL;
	if ((pfile = fopen (MiPfile, "rb")) == (FILE *)0) {
		printf("\nUnable to open file.\n");
		return 1;
	} 

	double tmp;
	long ii,jj,kk,r_1,c_1;

	// get number of symbols
	int k_sym; // number of symbols
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	k_sym = (int)round(tmp);

	// because this routine requires a binary alphabet, exit if it isn't
	if (k_sym != 2) {
		cout << "Requires a binary (2 symbol) alphabet." << endl;
		return 1;
	}

	// allocate memory for the di_vals shift amount array (d_1,...,d_k), and retrieve from file
	long *di_vals; // 1D array of integer shift amounts, (d_1,...,d_k)
	if (!(di_vals = new long[k_sym])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<k_sym; ii++) {
		if (fread(&tmp, sizeof(double), 1, pfile) < 1){
			printf("\nUnable to read from opened file.\n");
			fclose(pfile);
			return 1;
		}
		di_vals[ii] = (long)round(tmp);
	}

	// get simplex height (=matrix size for binary sequences, k=2)
	long hk; // simplex height, in units of 1/(k-1); for binary symbol set, hk is just the matrix size
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	hk = (long)round(tmp);

	// get delta
	double delta; // matrix regularity parameter
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	delta = tmp;

	// get regularity type (pad only, pad and splice, or splice only)
	int p_and_s; // type of regularity: 0 for pad only; 1 for pad and splice; 2 for splice only
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	p_and_s = (int)round(tmp);

	// get number of block rows in the M_i*P set
	long nbrows; // number of block rows in the M_i*P set
	if (fread(&tmp, sizeof(double), 1, pfile) < 1){
		printf("\nUnable to read from opened file.\n");
		fclose(pfile);
		return 1;
	}
	nbrows = (long)round(tmp);

	// allocate memory for mxsz_vals and read in matrix block size values
	long *mxsz_vals; // 1D array holding the block matrix sizes (nbrows total)
	if (!(mxsz_vals = new long[nbrows])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<nbrows; ii++) {
		if (fread(&tmp, sizeof(double), 1, pfile) < 1){
			printf("\nUnable to read from opened file.\n");
			fclose(pfile);
			return 1;
		}
		mxsz_vals[ii] = (long)round(tmp);
	}
	long mx_sz = mxsz_vals[0]; // matrix dimension--binary alphabet, so only one block row; this should also equal hk for the binary alphabet case



	//***
	// allocate memory for {M0*P,M_1*P,P} matrix set
	//***

	// there will be k_sym matrices for the nth block row (n from 1 to nbrows), each of the k_sym matrices being of dimension mxsz_vals[n-1]
	double **MiP; // 2D array for the matrix entry data; given block row i, symbol #j, matrix row r, column c, format is: MiP[i*k_sym+j][r*mxsz_vals[i]+c]
	if (!(MiP = new double*[nbrows*k_sym])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	for (ii=0; ii<nbrows; ii++) {
		for (jj=0; jj<k_sym; jj++) {
			if (!(MiP[ii*k_sym + jj] = new double[mxsz_vals[ii]*mxsz_vals[ii]])) { // assign pointer to appropriate allocated spot in memory
				cout << "Error: out of memory." << endl;
				exit(1);
			}
		}
	}

	// allocate memory for padding matrix P
	double *Pad;
	if (!(Pad = new double[mx_sz*mx_sz])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}



	//***
	// read in .bin file matrices
	//***

	// nb: this Mi*P reader is in a "universal format" (can handle multiple block rows); in this binary alphabet case, we'll only have 1 block row; the matrix indexing used below in the MC routine implicitly uses block row # = 0 (ie just MiP[jj][mx_sz*row+col], where jj is the symbol number and mx_sz = mxsz_vals[0] = hk)
	for (ii=0; ii<nbrows; ii++) {
		for (jj=0; jj<k_sym; jj++) {
			for (kk=0;kk<(mxsz_vals[ii]*mxsz_vals[ii]);kk++) {
	  			if (fread(&MiP[ii*k_sym+jj][kk], sizeof(double), 1, pfile) < 1) {
					printf("\nUnable to read from opened file.\n");
					return 1;
				}
			}
		}
	}

	// read in P (padding matrix); note in the splice-only case, this will just be the identity matrix
	for (jj=0;jj<(int)(mx_sz*mx_sz);jj++) {
		if (fread(&Pad[jj], sizeof(double), 1, pfile) < 1){
			printf("\nUnable to read from opened file.\n");
			return 1;
		}
	}

	fclose(pfile);



	//***
	// show header info
	//***

	cout << "Matrix data read:" << endl;
	cout << k_sym << " symbols; simplex height=" << hk << "; delta=" << delta;
	if (p_and_s==PAD_ONLY) {
		cout << "; pad only";
	} else if (p_and_s==PAD_AND_SPLICE) {
		cout << "; pad and splice";
	} else if (p_and_s==SPLICE_ONLY) {
		cout << "; splice only";
	} else {
		cout << "; p_and_s unknown";
	}
	cout << "; " << nbrows << " block rows." << endl;
	cout << "shift amounts: ";
	for (ii=0; ii<k_sym; ii++) {
		cout << di_vals[ii] << " ";
	}
	cout << endl;


	// DEBUG / for displaying / checking matrix loading
	/*
	ii = 0; // for which block row
	jj = 0; // for which symbol
	long mxsz_tmp = mxsz_vals[ii];
		for (r_1=0;r_1<mxsz_tmp;r_1++) { // row number
			for (c_1=0;c_1<mxsz_tmp;c_1++) { // column number
				cout << setprecision(8) << MiP[ii*k_sym+jj][mxsz_tmp*r_1+c_1] << " ";
			}
			cout << endl;
		}
		cout << endl;
	*/



	//***
	// allocate memory for Monte Carlo resampling
	//***

	// will need a 1-dimensional of length prod_len, to store the betas
	// will need 2 2-dimensional arrays, each of size num_sam by mx_sz, to store the (row) vectors that arise
		// from successive multiplications
	double l_t_nrm;
	double *alpha_vals = new double[num_sam];
	double *beta_vals = new double[prod_len];
	double *cu_alpha = new double[num_sam];
	double *c_N_vec = new double[mx_sz];
	double *c_N_tmp = new double[mx_sz];
	double *E_l; // will address row i (0 to num_sam-1), column j (0 to s-1) by E_l[i*s+j]
	if (!(E_l = new double[num_sam*mx_sz])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	double *F_l;
	if (!(F_l = new double[num_sam*mx_sz])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	double *MT_tot;
	if (!(MT_tot = new double[mx_sz*mx_sz])) {
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	double *log_h_vals_last = new double[num_sam];
	double *log_h_vals_curr = new double[num_sam];
	double *log_n_mult = new double[num_sam];
	double *log_nm_last = new double[num_sam];

	double *tmp_vec = new double[mx_sz];

	int l_1;
	int ix;
	double tmp_nrm;
	double norm_of_vec;
	double tmp_1;
	double c_N_log;
	double c_N_norm_log = 0.0;


	//srand(time(NULL));

	timeval t1;
	gettimeofday(&t1, NULL);
	srand(t1.tv_usec * t1.tv_sec); // this "fancy" timer is used in the case of rapid re-calls



	//***
	// initialize variables
	//***

	for (kk=0;kk<num_sam;kk++) {
		for (l_1=0;l_1<mx_sz;l_1++) {
			E_l[kk*mx_sz+l_1]=1; // ie we're keeping track of norm of matrix product by consolidating via operation
				// on column vector of all 1's
			F_l[kk*mx_sz+l_1]=0;
		}
		log_h_vals_last[kk] = 0;
		log_n_mult[kk] = 0;
		log_nm_last[kk] = 0;
	}

		// this puts (M_0+M_1)*P into MT_tot, and initializes norm computation vector c_N_vec
	for (r_1=0; r_1<mx_sz; r_1++) {
		for (c_1=0; c_1<mx_sz; c_1++) {
			MT_tot[mx_sz*r_1+c_1] = MiP[0][mx_sz*r_1+c_1] + MiP[1][mx_sz*r_1+c_1];
			//MT_tot[mx_sz*r_1+c_1] = MT_0[mx_sz*r_1+c_1] + MT_1[mx_sz*r_1+c_1];
		}
		c_N_vec[r_1] = 1.0;
	}



	//***
	// outer loop through symbols (bits); two inner loops through samples
	//***
		
	for (jj=0;jj<prod_len;jj++) {

		// compute normalization constant--reciprocal of ||P((M_0+M_1)P)^m||_{en}

			// takes past normalization vector, c_N_vec, and multiplies it (on left) by (M_0+M_1)Pad
		for (r_1=0; r_1<mx_sz; r_1++) {
			c_N_tmp[r_1] = 0.0;
			for (c_1=0; c_1<mx_sz; c_1++) {
				c_N_tmp[r_1] += MT_tot[mx_sz*r_1+c_1]*c_N_vec[c_1];
			}
		}

			// normalize and carry constant, to prevent overflows
		norm_of_vec = 0.0;
		for (r_1=0; r_1<mx_sz; r_1++) {
			norm_of_vec += c_N_tmp[r_1]*c_N_tmp[r_1];
		}
		norm_of_vec = sqrt(norm_of_vec);
		for (r_1=0; r_1<mx_sz; r_1++) {
			c_N_vec[r_1] = c_N_tmp[r_1]/norm_of_vec;
		}
		//cout << norm_of_vec << endl;
		c_N_norm_log += log(norm_of_vec);
		//cout << c_N_norm_log << endl;

			// c_N value
		tmp = 0.0;
		for (r_1=0; r_1<mx_sz; r_1++) {
			for (c_1=0; c_1<mx_sz; c_1++) {
				tmp += Pad[mx_sz*r_1+c_1]*c_N_vec[c_1]; // this is the entrywise norm of Pad((M_0+M_1)P)^jj (w/o cumulative norm of c_N_vec)
			}
		}
		c_N_log = log(tmp) + c_N_norm_log;


		for (kk=0;kk<num_sam;kk++) {

			// compute F_t = M_i P E_t
			if ((rand()%2)==1) { // M_0 P
				for (r_1=0; r_1<mx_sz; r_1++) {
					F_l[kk*mx_sz+r_1] = 0;
					for (c_1=0; c_1<mx_sz; c_1++) {
						F_l[kk*mx_sz+r_1] += MiP[0][mx_sz*r_1+c_1]*E_l[kk*mx_sz+c_1];
					}
				}
			} else { // M_1 P
				for (r_1=0; r_1<mx_sz; r_1++) {
					F_l[kk*mx_sz+r_1] = 0;
					for (c_1=0; c_1<mx_sz; c_1++) {
						F_l[kk*mx_sz+r_1] += MiP[1][mx_sz*r_1+c_1]*E_l[kk*mx_sz+c_1];
					}
				}
			}

			// compute entropy--in each (kkth sample) position of F_l[kk*mx_sz+l_1] (think of its kk-1th row),
				// is the vector representing the product M_i*P...M_i*P*[1...1]'--this needs a left
				// multiplication by Pad before the entrywise norm is computed
			for (r_1=0; r_1<mx_sz; r_1++) {
				tmp_vec[r_1] = 0.0;
				for (c_1=0; c_1<mx_sz; c_1++) {
					tmp_vec[r_1] += Pad[mx_sz*r_1+c_1]*F_l[kk*mx_sz+c_1];
				}
			}

			tmp_nrm = 0.0;
			for (l_1=0;l_1<mx_sz;l_1++) {
				tmp_nrm += tmp_vec[l_1]; // this gets the matrix s.o.e. norm, through the vector
			}
			l_t_nrm = log(tmp_nrm)+log_nm_last[kk]; // reconstitute the (log of the) true vector value
			log_h_vals_curr[kk] = l_t_nrm - c_N_log + log(c_N_log-l_t_nrm); // compute entropy
				//log(tmp_nrm/n) + log(log(n)-log(tmp_nrm)); i.e. log(-p log(p)) = log(p)+log(log(1/p))

			// (re)normalize the current F_l--this keeps values within reasonable bounds
			norm_of_vec = 0.0;
			for (l_1=0; l_1<mx_sz; l_1++) {
				norm_of_vec += F_l[kk*mx_sz+l_1]*F_l[kk*mx_sz+l_1];
			}
			norm_of_vec = sqrt(norm_of_vec);
			for (l_1=0; l_1<mx_sz; l_1++) {
				F_l[kk*mx_sz+l_1] = F_l[kk*mx_sz+l_1]/norm_of_vec; 
			}
			log_n_mult[kk] = log_nm_last[kk] + log(norm_of_vec);

		} // end 1st loop through samples
		/*
			to this point we have the vectors F_l (nb: these are not normalized--that only happens when the entropy is computed), and the log of the cumulative normalization factor, log_n_mult
			we have the log of the entropy, log_h_vals_curr (and the last step's log(entropies), log_h_vals_last)
		*/

		beta_vals[jj]=0;
		for (kk=0;kk<num_sam;kk++) {
			alpha_vals[kk] = exp(log_h_vals_curr[kk]-log_h_vals_last[kk]);
			beta_vals[jj]+=alpha_vals[kk];
			cu_alpha[kk] = beta_vals[jj];
		}

		for (kk=0;kk<num_sam;kk++)
			cu_alpha[kk] /= beta_vals[jj]; // if alpha_vals = 5,2,7, then beta_vals[jj]=14, cu_alpha = 5/14,1/2,1

		// resampling loop
		for (kk=0;kk<num_sam;kk++) {
			tmp = (rand()%1000000000)/1000000000.0; // in range [0,1)
			for (ix=0;ix<num_sam;ix++) {
				if (cu_alpha[ix]>tmp) {
					break;
				}
			}
			
			// copy F_l[ix] into E_l[kk]
			for (l_1=0;l_1<mx_sz;l_1++) {
				E_l[kk*mx_sz+l_1] = F_l[ix*mx_sz+l_1];
			}
			
			// copy log_n_mult[ix] into log_nm_last[kk]
			log_nm_last[kk] = log_n_mult[ix];

			// copy log_h_vals_curr[ix] into log_h_vals_last[kk]
			log_h_vals_last[kk] = log_h_vals_curr[ix];


		} // end resampling loop

		// keep track of not just main vectors on resampling, but the norm factors--this amounts to saving
			// log_n_mult values under the pre-resampled order, in log_nm_last values, in the post-resampled order
		// after resampling, save the entropy values calculated pre-resampling (log_h_vals_curr) in the correct
			// spots in the array log_h_vals_last

	} // end loop through bits



	//***
	// display entropy estimate
	//***

	tmp = 0.0;
	for (jj=0;jj<prod_len;jj++) {
		tmp += log(beta_vals[jj]);
	}
	tmp_1 = exp(prod_len*(log(2)-log(num_sam))+tmp);
	cout << "Entropy est:" << endl;
	printf("%.9f\n", tmp_1);



	//***
	// release memory and exit
	//***

	delete [] di_vals;
	di_vals = NULL;
	delete [] mxsz_vals;
	mxsz_vals = NULL;
	for (ii=0; ii<nbrows*k_sym; ii++) {
		delete [] MiP[ii];
	}
	delete [] MiP;
	MiP = NULL;
	delete [] Pad;
	Pad = NULL;
	delete [] c_N_vec;
	c_N_vec = NULL;
	delete [] c_N_tmp;
	c_N_tmp = NULL;
	delete [] alpha_vals;
	alpha_vals = NULL;
	delete [] beta_vals;
	beta_vals = NULL;
	delete [] cu_alpha;
	cu_alpha = NULL;
	delete [] E_l;
	E_l = NULL;
	delete [] F_l;
	F_l = NULL;
	delete [] MT_tot;
	MT_tot = NULL;
	delete [] log_h_vals_last;
	log_h_vals_last = NULL;
	delete [] log_h_vals_curr;
	log_h_vals_curr = NULL;
	delete [] log_n_mult;
	log_n_mult = NULL;
	delete [] log_nm_last;
	log_nm_last = NULL;
	delete [] tmp_vec;
	tmp_vec = NULL;


}



// for debugging
void cout_vec(double *vec, int s)
{
	int j;
	for (j=0; j<s; j++) {
		cout << vec[j] << " ";
	}
	cout << endl;
}




