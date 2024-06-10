# CL-SRDE algorithm for CEC 2024 constrained benchmark

A C++ implementation of the CL-SRDE algrithm (Constrained Linear population size reduction Success Rate-based adaptive Differential Evolution) for the Congress on Evolutionary Computation competition on single-objective numerical optimization (https://github.com/P-N-Suganthan/2024-CEC).

The algorithm code is in "CL-SRDE.cpp" file.

# Compilation and usage

Compilation is simple using gcc/g++:

g++ -std=c++11 -O3 CL-SRDE.cpp -o L-SRTDE.exe

or depending on hardware

g++ -std=c++11 -O3 -march=corei7-avx CL-SRDE.cpp -o CL-SRDE.exe

Please note that the compilation requires support of C++11 standard.

This will create CL-SRDE executable, available for running.

Data will be written to "CL-SRDE_F#_D#.txt", where F and D are the function number and problem dimention.
