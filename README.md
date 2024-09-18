# CL-SRDE algorithm for CEC 2024 constrained benchmark

A C++ implementation of the CL-SRDE algrithm (Constrained Linear population size reduction Success Rate-based adaptive Differential Evolution) for the Congress on Evolutionary Computation competition on single-objective numerical optimization (https://github.com/P-N-Suganthan/2024-CEC).

The algorithm code is in the "CL-SRDE.cpp" file.

# Compilation and usage

Compilation is simple using gcc/g++:

g++ -std=c++11 -O3 CL-SRDE.cpp -o L-SRTDE.exe

or depending on hardware

g++ -std=c++11 -O3 -march=corei7-avx CL-SRDE.cpp -o CL-SRDE.exe

Please note that the compilation requires support of C++11 standard.

This will create CL-SRDE executable, available for running. Extract inputData.zip before running algorithm.

# Results format

Data will be written to several files with different format. Here F and D are the function number and problem dimension:

1. CL-SRDE_f_F%d_D%d_.txt - stores the function values, 25x2000, 25 is the number of runs, 2000 steps of saving.

2. CL-SRDE_C_F%d_D%d.txt - stores constraint violation values, 25x18000, here for each of 2000 steps 9 values are saved. There are 9 because there are maximum 3 inequality constraints and 6 equality constraints in CEC 2017 benchmark. For example, for F1 there will be some value in the first place, and the other 8 will be 0: "-57252.2    0    0    0    0    0    0    0    0" and for F7 it will be "0    0    0    1865.03    -1865.03    0    0    0    0    0    0    0", as there are 2 equality constraints.

3. CL-SRDE_C2_F%d_D%d.txt - stores similar results to CL-SRDE_C_F%d_D%d.txt, but without zeros. Here the size of the file depends on the number of constraints, i.e. for F1 it will be 25x2000, as there is only one constraint, but for F7 it will be 25x4000.

4. CL-SRDE_F%d_D%d.txt - stores both the function values and total violations. The file is 50x2000, because there are pairs of lines. The first one is function value, and the second is the sum of all violations.

5. CL-SRDE_X%d_D%d.txt - stores the decision variables of the best solution. This file is 25x(2000*D), e.g. 25x60000 for D=30.

6. CL-SRDE_eps_F%d_D%d.txt - same as CL-SRDE_F%d_D%d.txt, but if total violation is below 0.0001, it is set to 0.
