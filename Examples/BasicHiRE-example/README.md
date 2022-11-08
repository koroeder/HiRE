1. Compiling HiRE

Requires: Fortran and C compilers and cmake. Here we use gfortran, gcc and cmake.

Set up a new directory 'build' in the pureHiRE directory, and enter this directory.

Run the command 'FC=gfortran CC=gcc cmake ../source' (replace gfortran and gcc with other compilers if used).
For a standard build, then run 'make'.

A successful compilation results in the executable HIRE.

For a debug version, before running make, use the command 'ccmake .' (requires cmake gui) and change the version from Release to Debug.

2. Running HIRE

All the input files are provided in the 'input' directory, copy the executable into the same directory as the files.
Run HIRE with the command './HIRE parameters.top scale_RNA.dat > output'.

To run with a different topology or parameter file, simply change the input names.
The name for the coordinate file is fixed as start.
