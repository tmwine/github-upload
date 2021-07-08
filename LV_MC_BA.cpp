/*

	Monte Carlo routine for estimating Lyapunov exponent and variance, for binary alphabets

	this expects, on the command line, <prod_len> <num_sam> <Mi*P data filename>
	<prod_len> is the number of matrices in each product chain
	<num_sam> is the number of product chains
	<Mi*P data file> is the name of the .bin file holding the relevant {M_1*P,M_2*P} matrix product pair

	this allows toggling between C++ <random> and C rand(), with the former a better randomizer but probably around 5-10% slower; to change the randomizer method:
		keyword search the code for "~RANDOMIZER~"
		un-comment lines with rand(), and comment out lines with <random> (or vice versa) to switch between the two types of randomization (rand() is faster, but is less uniform--for one, the remainder / wrap-around problem when wanting a uniform distribution in a range other than RAND_MAX--while <random> is closer to a uniform distribution, but around 10% slower)

	this reads in a standard format M_i*P matrix data file:
	<k_sym>, <di_vals>, <hk>, <delta>, <p_and_s>, <nbrows>, <mxsz_vals>, <MiP>
	k_sym = number of symbols
	di_vals = shift values (d_1,...,d_k)
	hk = simplex height (=matrix size in binary alphabet case)
	delta = pad/splice penalty parameter, for degree of regularity/randomness 
	p_and_s = 0 for pad only, 1 for pad and splice, 2 for splice only
	nbrows = number of block rows (=1 for binary alphabets)
	mxsz_vals = sizes of the matrices in respective block rows (array of elements nbrows long)
	MiP = matrix data; it's stored in order MiP_1,1,...,MiP_1,k; MiP_2,1,...,MiP_2,k; ...; ie the first block row, followed by the second block row, and so on
	nb: the individual MiP matrix entries are expected to be written in row-major order (that is, (r1,c1),(r1,c2),...,(r2,c1),(r2,c2),...,...)

 	nb: as of latest, the Octave routines calling this (via system call) require its first two lines of numeric output (strings / non-numeric lines notwithstanding) to be the LE est and the var est, respectively; eg the C++ routine output could be "LE est:\n1.235\nvar est:\n3.212", but not "LE est: 1.235\nvar est: 3.212"

	nb: this should compile on C++11 or later; eg $g++ -std=c++11 my_code.cpp -o my_code

*/


using namespace std; // this is generally bad practice

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <sys/time.h>
#include <iomanip> // for setprecision (at least)
#include <climits> // for checking bits in a byte
#include <random> // for optional C++ randomizer improvement over rand()

#define PAD_ONLY 0
#define PAD_AND_SPLICE 1
#define SPLICE_ONLY 2
#define RAND_RANGE 100000 // this is the desired resolution for the main Monte Carlo randomizer
#define RAND_ACC 10000.0 // this to check compatibility between RAND_RANGE and RAND_MAX (C library constant), for use with rand(); if RAND_MAX/RAND_RANGE is too small, the "remainder" error when modding rand()'s output by RAND_RANGE can unduly distort from uniform randomness
#define P_ERR 0.01 // this is a ~nice-to-have check for inconsistency in p_vals


void cout_vec(double *vec, int s); // for debugging


int main(int argc,char *argv[])
{

	if (argc != 4) {
		cout << "Problem with input arguments." << endl;
		exit(1);
	}

	long prod_len = atol(argv[1]);
	long num_sam = atol(argv[2]);
	char *myfile = argv[3];

	// check C++ here is using 8 bit bytes, and 8 bytes to a "double" (for Octave .bin file compatability; as of this writing)
	if (CHAR_BIT!=8) {
		cout << "Number of bits in a byte here is " << CHAR_BIT << ", not 8. This will likely produce problems with reading the Octave-generated M_i*P .bin data file, which renders doubles with 64 bits. Exiting." << endl;
		exit(1);
	} else if (sizeof(double)!=8) {
		cout << "C++ defines type double as " << sizeof(double) << " bytes. This will likely produce problems reading the Octave-generated M_i*P .bin data file, which renders doubles with 64 bits. Exiting." << endl;
	}



	//***
	// read in .bin file header
	//***

	// get prelim parameters for matrix set
	// the matrix data file has (in order) {matrices' size, number of symbols, p_1,...,p_k, M_1,...,M_k}

	FILE *pfile = NULL;
	if ((pfile = fopen (myfile, "rb")) == (FILE *)0) {
		printf("\nUnable to open file.\n");
		return 1;
	} 

	double tmp;
	long ii,jj,kk,rr,cc;
	long rand_mod = RAND_RANGE; // this assignment for easier toggling between the C rand()  version and the C++ <random> version

	// for ~RANDOMIZER~ RAND()
	//if (RAND_RANGE > RAND_MAX/RAND_ACC) { // check that RAND_RANGE isn't too large relative to RAND_MAX (remainder error)
	//	rand_mod = RAND_MAX;
	//	cout << "RAND_RANGE and RAND_MAX mismatch for Monte Carlo routine; defaulting to RAND_MAX" << endl;
	//}

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

	// allocate memory for MiP arrays, to hold matrices; there will be k_sym matrices for the nth block row (n from 1 to nbrows), each of the k_sym matrices being of dimension mxsz_vals[n-1]
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


	//***
	// read in .bin file matrices
	//***

	// nb: this reader is in a "universal format" (can handle multiple block rows); in this binary alphabet case, we'll only have 1 block row; the matrix indexing used below in the MC routine implicitly uses block row # = 0 (ie just MiP[jj][mx_sz*row+col], where jj is the symbol number and mx_sz = mxsz_vals[0] = hk)
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
	cout << "; pl=" << prod_len << "; ns=" << num_sam;
	cout << endl;

	cout << "block row sizes: ";
	for (ii=0; ii<nbrows; ii++) {
		cout << mxsz_vals[ii] << " ";
	}
	cout << endl;


	// DEBUG / for displaying / checking matrix loading
	/*
	ii = 0; // for which block row
	jj = 0; // for which symbol
	long mxsz_tmp = mxsz_vals[ii];
		for (rr=0;rr<mxsz_tmp;rr++) { // row number
			for (cc=0;cc<mxsz_tmp;cc++) { // column number
				cout << setprecision(8) << MiP[ii*k_sym+jj][mxsz_tmp*rr+cc] << " ";
			}
			cout << endl;
		}
		cout << endl;
	*/
	//


	//***
	// create pdf for Bernoulli matrix selection
	//***

		// create p_vals vector, and fill it
	double *p_vals;
	if (!(p_vals = new double[k_sym])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	double tmp_1 = 0.0;
	for (ii=0; ii<k_sym; ii++) {
		tmp_1 += 1/(float)di_vals[ii];
	}
	for (ii=0; ii<k_sym; ii++) {
		p_vals[ii] = (1/(float)di_vals[ii])/tmp_1;
	}

		// create the cdf of the p_vals distribution (for Bernoulli Monte Carlo matrix selection)
	double *pv_cdf;
	if (!(pv_cdf = new double[k_sym])) { // assign pointer to appropriate allocated spot in memory
		cout << "Error: out of memory." << endl;
		exit(1);
	}
	pv_cdf[0] = p_vals[0];
	for (ii=1; ii<k_sym; ii++) {
		pv_cdf[ii] = pv_cdf[ii-1]+p_vals[ii];
	}
	if (fabs(pv_cdf[k_sym-1]-1)>P_ERR) {
		cout << "Error: p_vals in file do not sum closely to 1." << endl;
		return 1;
	}


	//***
	// start Monte Carlo routine
	//***

	double acc_le = 0.0;
	double acc_le_tot = 0.0;
	double acc_les = 0.0;

	long mx_sz = mxsz_vals[0]; // binary alphabet, so only one block row; this should also equal hk for the binary alphabet case
	double *v_0 = new double[mx_sz];
	double *v_tmp = new double[mx_sz];

	int r_1,c_1;
	double norm_v_0;
	double rand_val;

	// ~RANDOMIZER~ RAND()
	//srand(time(NULL));

	// ~RANDOMIZER~ <random>
	random_device rd;  // for obtaining a random seed
	mt19937 gen(rd()); // standard mersenne_twister_engine seeded with rd()
	uniform_int_distribution<> dis(0,rand_mod);


	for (kk=0;kk<num_sam;kk++) {

		// make normalized random direction vector
		for (r_1 = 0; r_1<mx_sz; r_1++) {
			//v_0[r_1] = (rand()%(rand_mod+1))/(double)rand_mod; // ~RANDOMIZER~ RAND()
			v_0[r_1] = dis(gen)/(double)rand_mod; // ~RANDOMIZER~ <random>
		}
		norm_v_0 = 0.0;
		for (r_1 = 0; r_1<mx_sz; r_1++) {
			norm_v_0 += v_0[r_1]*v_0[r_1];
		}
		norm_v_0 = sqrt(norm_v_0);
		for (r_1 = 0; r_1<mx_sz; r_1++) {
			v_0[r_1] = v_0[r_1]/norm_v_0;
		}

		// DEBUG
		//cout_vec(v_0,mx_sz);

		acc_le = 0.0;

		for (ii=0; ii<prod_len; ii++) {
			//rand_val = (rand()%(rand_mod+1))/(double)rand_mod; // ~RANDOMIZER~ RAND()
			rand_val = dis(gen)/(double)rand_mod; // ~RANDOMIZER~ <random>

			// generalized Bernoulli Monte Carlo matrix selection
			for (jj=0; jj<k_sym; jj++) {
				if (rand_val <= pv_cdf[jj])
					break;
			}

			// uniform Monte Carlo matrix selection
			//jj = rand()%k_sym;

			// perform the matrix multiplication (in the chain M_i P M_i P ... M_i P)
			for (r_1=0; r_1<mx_sz; r_1++) {
				v_tmp[r_1] = 0;
				for (c_1=0; c_1<mx_sz; c_1++) {
					v_tmp[r_1] += MiP[jj][mx_sz*r_1+c_1]*v_0[c_1];
				}
			}

			norm_v_0 = 0.0;
			for (r_1=0; r_1<mx_sz; r_1++) {
				norm_v_0 += v_tmp[r_1]*v_tmp[r_1];
			}
			norm_v_0 = sqrt(norm_v_0);
			for (r_1=0; r_1<mx_sz; r_1++) {
				v_0[r_1] = v_tmp[r_1]/norm_v_0;
			}
			tmp = log(norm_v_0);
			acc_le += tmp;
		} // end loop over number of matrices in this product (ii)

		tmp_1 = acc_le/prod_len; // this is this run of prod_len matrices' Lyapunov exponent estimate
		acc_le_tot += acc_le;
		acc_les += tmp_1*tmp_1; // ie this takes the Lyapunov exponent estimate from this run of prod_len matrices, and squares it

	} // end loop over number of matrix product samples (kk)



	//***
	// display results
	//***

	cout << "Lyapunov exponent est:" << endl;
	printf("%.9f\n", acc_le_tot/(((double)prod_len)*num_sam));
	cout << "Lyapunov variance est:" << endl;
	printf("%.9f\n", prod_len*(acc_les/num_sam-pow(acc_le_tot/(((double)prod_len)*num_sam),2)));



	//***
	// release memory and exit
	//***

	delete [] v_0;
	v_0 = NULL;
	delete [] v_tmp;
	v_tmp = NULL;
	delete [] p_vals;
	p_vals = NULL;
	delete [] pv_cdf;
	pv_cdf = NULL;
	delete [] di_vals;
	di_vals = NULL;
	delete [] mxsz_vals;
	mxsz_vals = NULL;
	for (ii=0; ii<nbrows*k_sym; ii++) {
		delete [] MiP[ii];
	}
	delete [] MiP;
	MiP = NULL;

	return 0;
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




