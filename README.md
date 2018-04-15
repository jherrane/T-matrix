The code  ("T_matrix") computes a T-matrix of an arbitrarily shaped inhomogeneous dielectric object by the electric current volume integral equation method (JVIE).

Requirements:
FFTW3
HDF5
LAPACK


Usage:
./T_matrix -parameter1 value1 -parameter2 value2 ...

-parameters default value
-k 1.0            wavenumber 
-mesh mesh.h5     input mesh file
-elem_ka 10.0     size parameter of the circumscribing sphere
-T_out T.h5       output T-matrix file


Another code ("multi_T") compiles an input file T_multi_ic.h5 for the R2T2 solver or T_multi_tot.h5 for the FaSTMM solver from the T-matrices computed with the JVIE code.

Usage ./multi_T -parameter1 value1 -parameter2 value2 ...

-parameters default value
-elem_ka 10.0                size parameter of circumscribing sphere
-N_Tin 1                     number of input T-matrices
-T_in T                      prefix for input T-matrices (T_1.h5, T_2.h5, ..., T_{N_Tin})
-T_multi T_multi_tot.h5      output file1 (input for FaSTMM)
-T_multi_ic T_multi_ic.h5    output file2 (input for R2T2)