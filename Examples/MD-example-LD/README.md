1. Compiling the MD code

Requires: Fortran and C compilers and cmake. Here we use gfortran, gcc and cmake.

Set up a new directory 'build' in the MD_HiRE directory, and enter this directory.

Run the command 'FC=gfortran CC=gcc cmake ../source' (replace gfortran and gcc with other compilers if used).
For a standard build, then run 'make -j8' (using 8 cores in parallel, a plain make also works).

A successful compilation results in the executable HIREMD.

For a debug version, before running make, use the command 'ccmake .' (requires cmake gui) and change the version from Release to Debug.

2. Running the MD code.

All the input files are provided in the 'input' directory - just run the executable HIREMD in any directopry with those files.

Options:

 - to change the timestep alter the line TIMESTEP 0.01 to a different value
 - adjust the Langevin dynamics settings with changes to GAMMA 0.1
 - for more MD steps, change the MDSTEPS value
 - for a different temperature, change the line TEMPERATURE 0.616 (the units are kcal/mol!)

