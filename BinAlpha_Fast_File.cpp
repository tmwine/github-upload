/*

	fast matrix multiplication, for binary alphabets
	text file string input (individual symbols in range [1,...,k_sym])

	this expects, on the command line, <Mi*P data filename> <symbol string filename>
	<Mi*P data filename> is the appropriate .bin file containing the relevant Mi*P product pair
	<symbol string filename> should be a text file, containing only lines of numbers (no delimiting--no commas or spaces), each in the range [0,1] (eg "011010101001"); nb: this expects the values in symbol_string to be in the valid range--it does not check this beforehand--if there is an error, it will give an error message
	nb: as of latest, this assumes a consecutive mapping between chars '0','1',...'9' and their representative (byte) values--this to quickly convert a char decimal digit into its value (through char-'0')
	nb: the value of LINES_PER_BLOCK in the header area sets how many lines to read successively from the file before performing a matrix multiplication (processes in batches)--this may be adjusted as necessary

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

	the value of LINES_PER_BLOCK in ~header sets how many lines to read successively from the file before performing a matrix multiplication (processes in batches)

	nb: this internally changes the multiplication order from the command line string input version, which was S_m*...*S_1*[1,...,1]', where the multiplication necessarily went from right to left (which required reading the input string backwards, to form the product; reading backwards is not that practical when reading a text file, especially a large one); this version multiplies from left to right, via [1,...,1]*S_m*...*S_1 (ie v_0 has been changed from a column vector to a row vector, and where by assoc. the sum of elements in the resulting vector still gives the entrywise norm)

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

	if (argc != 3) {
		cout << "Problem with input arguments." << endl;
		exit(1);
	}

	char *MiPfile = argv[1];
	char *symbol_string_file = argv[2];

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

	// get prelim parameters for matrix set
	// the matrix data file has (in order) {matrices' size, number of symbols, p_1,...,p_k, M_1,...,M_k}

	FILE *pfile = NULL;
	if ((pfile = fopen (MiPfile, "rb")) == (FILE *)0) {
		printf("\nUnable to open file.\n");
		return 1;
	} 

	double tmp;
	long ii,jj,kk,rr,cc;

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


	//***
	// initialize and set up file reading, and allocate memory
	//***

	long r_1,c_1;
	double norm_v_0;
	double log_vec_norm = 0.0;
	int sym_val;
	bool bad_char = false;

	long mx_sz = mxsz_vals[0]; // binary alphabet, so only one block row; this should also equal hk for the binary alphabet case
	double *v_0 = new double[mx_sz]; // since this routine multiplies "left to right," and we want v_0*M_i*P*...*M_i*P, where v_0 is a row vector, v_0 here is thought of as a row vector
	double *v_tmp = new double[mx_sz]; // temp (row) vector, for intermediary calculations
	for (c_1 = 0; c_1<mx_sz; c_1++) {
		v_0[c_1] = 1.0; // sets up the left-hand-side row vector of 1's, to start the multiplication (left to right)
	}

	ifstream string_file(symbol_string_file);
	bool file_done = false;
	string file_line;
	string line_block;




	//***
	// matrix multiplication routine
	//***

	// nb: this changes the multiplication order from BinAlpha_Fast_CL (command line symbol string input version)--BinAlpha_Fast_CL was S_m*...*S_1*[1,...,1]', where the multiplication necessarily went from right to left (which required reading the input string backwards, to form the product; reading backwards is not that practical when reading a (possibly large) text file); this version multiplies from left to right, via [1,...,1]S_m*...*S_1 (ie v_0 has been changed from a column vector to a row vector, and where by assoc. the sum of elements in the vector still gives the entrywise norm)

	while (!file_done) {

		line_block = "";
		for (jj=0;jj<LINES_PER_BLOCK;jj++) {
			if (getline(string_file, file_line)) {
				line_block += file_line;
			} else {
				file_done = true;
				break;
			}
		}

		// loop through symbols in line_block string
		// this keeps track of the norm of v_0*MiP...M_iPv_0, where v_0 is initially [1...1]
		for (kk=0;kk<line_block.length();kk++) {

			sym_val = line_block[kk]-'0';

			if (sym_val<0 || sym_val>=k_sym) {
				bad_char = true;
				break;
			}

			// perform the matrix multiplication (in the chain M_i P M_i P ... M_i P [1,...,1]')
			for (c_1=0; c_1<mx_sz; c_1++) {
				v_tmp[c_1] = 0.0;
				for (r_1=0; r_1<mx_sz; r_1++) {
					v_tmp[c_1] += v_0[r_1]*MiP[sym_val][mx_sz*r_1+c_1];
				}
			}

			// DEBUG
			//cout_vec(v_0,mx_sz);

			norm_v_0 = 0.0;
			for (c_1=0; c_1<mx_sz; c_1++) {
				norm_v_0 += v_tmp[c_1]*v_tmp[c_1];
			}
			norm_v_0 = sqrt(norm_v_0);
			for (c_1=0; c_1<mx_sz; c_1++) {
				v_0[c_1] = v_tmp[c_1]/norm_v_0;
			}
			log_vec_norm += log(norm_v_0);

		} // end loop over symbols in symbol string


	} // end while loop over all the lines in the symbol string file



	//***
	// output results
	//***

	// at this point, v_0 stores the normed result of M_i P M_i P ... M_i P [1,...,1]' (ie a vector, or array of length mx_sz), with log of the norm in log_vec_norm

	double log_en_nrm = 0.0;

	if (bad_char) {
		cout << "symbol out of range in input string; multiplication failed" << endl;
	} else {
		// reconstitute the log of the entrywise norm
		double tmp_sum = 0.0;
		for (c_1=0; c_1<mx_sz; c_1++) {
			tmp_sum += v_0[c_1];
		}
		log_en_nrm = log_vec_norm + log(tmp_sum);
		cout << "log of entrywise norm:" << endl;
		printf("%.14f\n",log_en_nrm);
	}


	//***
	// release memory and exit
	//***

	delete [] v_0;
	v_0 = NULL;
	delete [] v_tmp;
	v_tmp = NULL;
	delete [] di_vals;
	di_vals = NULL;
	delete [] mxsz_vals;
	mxsz_vals = NULL;
	for (ii=0; ii<nbrows*k_sym; ii++) {
		delete [] MiP[ii];
	}
	delete [] MiP;
	MiP = NULL;



	if (!bad_char) {
		return 0;
	} else {
		return 1;
	}


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




