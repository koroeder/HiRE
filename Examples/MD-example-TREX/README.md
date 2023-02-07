1. Compiling the MD code with MPI

Requires: Fortran and C MPI compilers, openmpi (or another mpi version) and cmake. Here we use mpif90 (gfortran), mpicc (gcc) and cmake.

Set up a new directory 'build/MPI-HiRE' in the MD_HiRE directory, and enter this directory.

Run the command 'FC=mpif90 CC=mpicc cmake ../source -DCOMPILER_SWITCH=gfortran -DWITH_MPI=ON' (replace with other compilers if used).
The DWITH_MPI flag can also be set with ccmake using a GUI, but the DCOMPILER_SWITCH flag is required to generate all files required by make.
For a standard build, then run 'make -j8' (using 8 cores in parallel, a plain make also works).

A successful compilation results in the executable HIREMD.

For a debug version, before running make, use the command 'ccmake .' (requires cmake gui) and change the version from Release to Debug.

2. Running the MD code.

All the input files are provided in the 'input' directory.
To run the MPI version of HIREMD, run 'mpirun -np 4 ${PATHTO}/HIREMD'.
Replace ${PATHTO} with your path to the hire executable. The argument after -np is the number of process (cores).
This will need to match the number of replicas in mddata!

Options:
 - number of replicas and temperature range (REXMD line)
 - all other MD options available in plain MD simulations
 - change "T" to "H" for Hamiltonian-REX MD
