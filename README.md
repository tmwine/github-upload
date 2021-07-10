# github-upload

This repository contains the routines for assigning a regularity z-score to a symbolic sequence, as well as determining the entropy of associated distributions in the binary case (from the paper [A Matrix-Based Regularity Measure for Symbolic Sequences](https://osf.io/vpg8h)). Most are written in Octave. A few are in C++, for the sake of speed.

First are some examples. The routines and main variable types are explained in more detail afterward.

### Some examples using the routines

These examples assume the .cpp files have been compiled under their respective filenames.

  #### z-score under a binary alphabet

Assign a z-score to the sequence 00100010101000011111 under p=1/2, matrix dim=2, delta=0.5, pad-and-splice.

In Octave, run 
```
>>BinAlpha_Mx_Gen(2,1,1,0.5,1)
```
and save when prompted (say as MiP_file.bin).

If the Lyapunov exponent and variance have not yet been estimated,
- from the command line, run

    ```
    $./LV_MC_BA prod_len num_sam MiP_file.bin
    ```
    (let lambda be the Lyapunov exponent, and varn the variance; lambda~-0.417725, varn~0.0255641)
- (alternatively, hard code the name of the binary Monte Carlo .cpp file into BinAlpha_Le_Var_MC.m; then in Octave, run
    ``` 
    >>BinAlpha_Le_Var_MC(prod_len,num_sam,MiP_file.bin)
    ```
    )
- (nb: parameter values for robust results from brute force Monte Carlo with binary alphabets tended to be of order prod_len=1000000, num_sam=100000 (provided matrix size is not too large, say < 50))
- as an optional step, sample vs theoretical bell curves can be plotted, and a correction factor can be estimated by running in Octave

    ```
    >>Bell_Plot(prod_len,num_sam,lambda,varn,MiP_file.bin)
    ```
    where prod_len is approximately the product length under test (here, 20), and num_sam is sufficiently large to get a smooth sample curve (in practice, say of order 1e4); the correction factor (call it cf) becomes more important when the sequence length is short, of order 10*(matrix dim) (as a rule of thumb, the sequence under test should not be shorter); in this example, cf is estimated to be 0.7

From the command line, run
```
$./BinAlpha_Fast_CL MiP_file.bin 00100010101000011111
```
This gives x=log||.||=-7.93455. Then compute the z-score: z-score = (x-(prod_len\*lambda+cf))/sqrt(prod_len\*varn) = -0.38764.

  #### entropy of a given distribution

Estimate the entropy of binary sequences of length 50 under the distribution with parameters p=3/5, matrix dim=5, delta=0.25, pad-only.

In Octave, run 
```
>>BinAlpha_Mx_Gen(5,2,3,0.25,0)
```
and save when prompted (say as MiP_file.bin).

From the command line, run 
```
$./BinAlpha_Entropy prod_len num_sam MiP_file.bin
```
(nb: prod_len must be kept short for this to execute in reasonable timespans (on typical 2020's-era processors, say < 200), and num_sam of order 1e5).)

Note, from the results of the paper, if h_est is the result of the resampled Monte Carlo routine, then the average entropy, h_est/prod_len, should be approximately constant. Also, for alphabets of > 2 symbols there is not yet an entropy routine; this would just involve incorporating block row matrix multiplication in the BinAlpha_Entropy routine
 

  #### z-score under an alphabet of 3 symbols

Assign a z-score to the 300-long 3-symbol sequence
```
231112112222121122213111212111121212221221213111223312311112211121232121233323131221222333312312112132311313111121122122211123221311311231112212312132233122232213122211121122231111212211212311111231112232222221121122312211211211113112322112211122212121123122233121121131131232311212332333121131222221
```
under shift amounts [2,2,5] (corresponding to probabilities P(1)=P(2)=5/12, P(3)=1/6), at simplex height 16 (minimal), delta=0.01, pad-and-splice.

In Octave, run 
```
>>MultAlpha_Mx_Gen(3,[2,2,5],16,0.01,1)
```
and save the Si forms (say as Si_file) and the Mi*P forms (say as MiP_file.bin) when prompted. This produces 5 block rows: 2 of dimension 18, 2 of dimension 19, and 1 of dimension 20.

To compute the Lyapunov exponent and variance, there are two main options:
- option 1 (low level): if nbrows is the number of block rows, run from the command line
    ```
    $./LV_MC_MA prod_len num_sam MiP_file.bin block_row_#
	``` 
    for each block row. Note, if there are symmetries (for instance if some of the d_i are equal to each other), this can be unecessarily time consuming because LE/var values for isomorphic block rows are necessarily the same.
- option 2 (checking for symmetries first): in Octave, run
    ```
    >>MultAlpha_Iso_Fam("Si_file")
    ```
    and save the output block row family file (say as fam_file); then from Octave run
    ```
    >>MultAlpha_Le_Var_MC(prod_len,num_sam,"MiP_file.bin","fam_file")
    ```
    This has the advantage of MultAlpha_Iso_Fam having checked for block row isomorphisms, so MultAlpha_Le_Var_MC's internal calls to LV_MC_MA will be run only on necessary block rows (mutually non-isomorphic). This shows the 2 block rows of dimension 18 are isomorphic, and so are the 2 block rows of dimension 19, making 3 block row families total.

For the LE/var estimates: #1, LE: -2.033874, var: 3.633530; #2, LE: -1.914878, var: 3.926234; #3, LE: -1.877799, var: 4.010124.

The Lyapunov exponents from the different families show minimal separation magnitude of ~0.04. Z-score calculation will use the maximal Lyapunov exponent, of -1.877799. For z-scores to be accurate, the length of the sequence under test should be at least as long as 10 times the matrix dimension--so > 200 symbols. The other consideration is, ideally, sufficient separation from the other block rows' bell curves. At a length m=300, the mean of the largest and next-to-largest block row families will be around -563 and -574, respectively. With a standard deviation of ~2 for each curve, this should allow for sufficient separation.

As an optional step, the LE/var calculations can be checked, and a correction factor can be estimated by running in Octave 
```
>>Bell_Plot(prod_len,num_sam,lambda,varn,"MiP_file.bin")
```
where lambda and varn correspond to the block row with the largest Lyapunov exponent (in this case, that of block row family #3), and prod_len is approximately the product length under test, and num_sam is sufficiently large to get a relatively smooth sample curve; the correction factor (call it cf) becomes more important when the sequence length is short; in this example, cf is estimated to be 9.0.

As another optional step, the block row Si template matrices' graphs can be displayed. From Octave, run
```
MultAlpha_Graph_Disp("Si_file",block_row_#,"output_file")
```
where output_file is used as an intermediary Graphviz file. Then process with Graphviz's Kneato--from the command line, run
```
$dot -Kneato -Gsplines=true -Gsep=.3 -Tps output_file -o graphout.ps
```
Here is the cylinder graph corresponding to the first block row:
![image](https://github.com/tmwine/github-upload/blob/master/pictures/graph_2_2_5_16.jpg)


To complete the z-score estimation, from the command line, run
```
$./MultAlpha_Fast_CL MiP_file.bin 23111...
```
This gives x = log||.|| = -561.8678. Then compute the z-score: z-score = (x-(prod_len\*lambda+cf))/sqrt(prod_len\*varn) = -0.2170443138.



### Variable types
   
 
- data_header--an Octave structure describing the sequence type (and regularity template information); its fields are:
  - .k_sym (scalar): number of symbols
  - .di_vals (vector): jth component is the geometric shift amount for symbol j
  - .hk (scalar): height of simplex (for binary sequences, this is the number of positions available for the one-dimensional walk)
  - .delta (scalar): randomness/regularity parameter, in [0,1]: 0 for pure regularity, 1 for pure (generalized Bernoulli) randomness
  - .p_and_s (scalar): the type of sequence operations permitted to match the regularity templates, in {0,1,2}: 0 for pad-only; 1 for pad-and-splice; 2 for splice-only
  - .br_dims (vector): block row matrix dimensions (for binary alphabets, this vector has one element, its value matching the value in data_header.hk)
      
- nbrows (scalar): the number of block rows in the matrix set
    
- Si_mx_cell--an Octave cell array of size [nbrows,k_sym]; the ith row contains the k_sym binary-valued nilpotent template matrices corresponding to the ith block row of the {S<sub>1</sub>,...,S<sub>k_sym</sub>} shift matrix set for the sequence in question
    
- graph--an Octave structure; its fields hold detailed information on the graph corresponding to a specific block row of an S_i shift matrix set; its fields are
	- .cyl_graph_matrix (matrix):
      - a temporal cylindrical graph matrix, the result of modding the graph's time steps by min_orbit; the resulting graph has exactly min_orbit time steps, the nodes at time step min_orbit always connecting to one or more nodes at time step 1
      - each row corresponds to a single node, its entries as follows: [time step, this node # (at this time step), exit edge symbol, edge-to-exit-node original (tot_mix) #, this node original #, time step of exit node (either the next time step, or a loopback), node # of exit node (at that time step)]
	- .binary_edge_matrix (vector of strings, columnar):
      - used for isomorphism checking
      - there are max_node_count\*min_orbit rows total, each row a binary string of 2\*k_sym characters; the ith time step in the graph corresponds to rows max_node_count\*(i-1)+1 to max_node_count\*i; time step i, node j (at that time step) corresponds to row max_node_count\*(i-1)+j (some rows are not used--they are padded with strings of all "0"'s)
      - the row string binary codes the edges entering the node, in the first k_sym characters, and the edges leaving the node, in the last k_sym characters--the jth character (left to right) in the string half correponds to the presence ('1') or absence ('0') of the jth symbol (in the edges entering, or leaving that node)--eg "010110" with k_sym=3 indicates an edge with symbol #2 enters the node, and edges with symbol #1 and #2 leaves the node
	- .min_orbit (scalar): the (time step) length used to make the graph cylindrical (ie time mod min_orbit); this value should always match the highest time step achieved in the respective temporal node matrix (eg last value in 1st column of temporal node matrix)
	- .max_node_count (scalar): stores the maximum number of nodes over all time steps in the post-min_orbit-modded temporal graph form; this value is important for parsing the respective binary edge matrix; the following equality should always hold: max_node_count\*min_orbit = # of rows in respective binary edge matrix
	- .node_count_vec (vector): has min_orbit elements, the ith entry storing the node count at time step i
	- .node_numbering_vec (vector, columnar): holds the original node numbers; it matches, row for row, the nodes referenced in the corresponding binary edge matrix; note: the binary edge graph will have the binary edge in/out strings ordered in (decreasing) lexicographical order, and the node_numbering_vec will match that ordering
	- .min_multiplicity_node (scalar): saves the node number of a node with binary edge code of minimal multiplicity for the graph
    
- graph_array: an Octave array of structures; the ith entry is an Octave graph structure corresponding to the ith block row of the S_i matrix set
    
- br_iso_families:
  - an Octave cell array of cell arrays; there is one outer cell array for each matrix dimension (in increasing order of dimension), and one inner cell array for each family (of block row indices), at that matrix dimension
  - eg for 5 block rows, if block rows 1,3, and 5 are of dimension 20, and block rows 2 and 4 are dimension 15, with 1 and 3 isomorphic (and the rest not), then br_iso_families{1} corresponds to block rows of dimension 15, and br_iso_families{2} corredsponds to block rows of dimension 20; br_iso_families{1} holds the cell vector of length 2, the first cell a single element vector, [2], and the second cell also a vector with one element, [4]; br_iso_families{2} holds a cell vector of length 2, the first cell a vector of length 2--[1,3], and the second cell a vector of length 1--[5]
    
- fam_LE_var_mx (matrix): used to record the results of the Monte Carlo Lyapunov exponent and variance estimates; each row corresponds to an indexed isomorphism family # (the order corresponding to order of appearance in br_iso_families); each row has elements [LE, var, product length, number of samples]


### Routines

  #### binary alphabets:
  - BinAlpha_Mx_Gen.m: generates the {M<sub>1</sub>\*P,M<sub>2</sub>\*P,P} matrix set; gives option to save result in a .bin file
  - BinAlpha_Le_Var_MC.m: manages the Monte Carlo estimates for Lyapunov exponent and variance; drives LV_MC_BA.cpp's executable via system call
      - LV_MC_BA.cpp: Monte Carlo estimates for Lyapunov exponent and variance; receives an Mi\*P .bin filename
  - BinAlpha_Fast_CL.cpp: fast matrix multiplier routine; receives an Mi\*P .bin filename; enter binary symbol sequence on the command line; returns the log of the entrywise norm of the product
  - BinAlpha_Fast_File.cpp: fast matrix multiplier routine; receives an Mi\*P .bin filename; reads binary symbol sequence from text file; returns the log of the entrywise norm of the product
  - BinAlpha_Entropy.cpp: entropy estimator for binary distributions; receives an Mi\*P .bin filename and returns the entropy estimate associated with sequences of a fixed length


  #### alphabets with > 2 symbols:
  - MultAlpha_Mx_Gen.m: generates Si_mx_cells, graph_array, and Mi\*P matrix set; uses MultAlpha_Pare and MultAlpha_Routes; option to save results to standard files
    - MultAlpha_Pare.m: used to pare rows and columns of zeros from initial simplex shift matrices
    - MultAlpha_Routes.m: finds all cyclical routes in the pared matrix set; each set of connected cycles is portioned into a block row; uses MultAlpha_Graph_Gen
        - MultAlpha_Graph_Gen.m: generates the graph_array array of Octave graph structures
  - MultAlpha_Iso_Fam.m: given a .m file containing Si_mx_cells and graph_array data, finds all isomorphisms between block rows of the same matrix dimension, sorting them into families in a br_iso_families cell array; uses Iso_Perm_Layer and Iso_Full_Layer; option to save results
    - Iso_Perm_Layer.m: initial check for potential isomorphism between two graphs, a test graph compared against a base graph; uses Iso_LH_Graph_Chk
    - Iso_Full_Layer.m: completion of check of potential isomorphism between two graphs, test and base pairs that passed the checks at Iso_Perm_Layer; uses Iso_LH_Graph_Chk
        - Iso_LH_Graph_Chk: creates a canonical node matrix for a graph, dependent on the starting node; uses a left-hand-turn rule
  - MultAlpha_Le_Var_MC.m: manages the Monte Carlo estimates for Lyapunov exponent and variance; drives LV_MC_MA.cpp's executable via system call
      - LV_MC_MA.cpp: Monte Carlo estimates for Lyapunov exponent and variance; receives an Mi\*P .bin filename
  - MultAlpha_Fast_CL.cpp: fast matrix multiplier routine; receives an Mi\*P .bin filename; enter binary symbol sequence on the command line; returns the log of the entrywise norm of the product
  - MultAlpha_Fast_File.cpp: fast matrix multiplier routine; receives an Mi\*P .bin filename; reads binary symbol sequence from text file; returns the log of the entrywise norm of the product
  - MultAlpha_Graph_Disp.m: for alphabets of > 2 symbols, given an .m file containing Si_mx_cells and graph_array data, creates a .dot (.gv) file that can be rendered in eg GraphViz to display a given block row's graph

  #### for any alphabet:
  - Slow_Mult.m: a slower (Octave) routine for multiplying the matrices corresponding to an input symbol sequence; returns the log of the entrywise norm of the product
  - Bell_Plot.m: for a given sample size, will compute (generalized Bernoulli) random matrix products and plot their histogram; the plot is compared against the theoretical Gaussian, from the Lyapunov exponent and variance, in accordance with the central limit theorem; will also estimate a correction factor--for obtaining more accurate z-scores for relatively short symbolic sequences

### License

This project is licensed under the terms of the MIT license.


