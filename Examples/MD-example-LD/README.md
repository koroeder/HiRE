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

 - to change the timestep alter the line TIMESTEP 0.0005 to a different value (the units are ps!)
 - adjust the Langevin dynamics settings with changes to GAMMA 2.0 (the units are ps^-1!)
 - to give TIMESTEP and GAMMA in the internal time unit instead (sqrt(amu*A^2/(kcal/mol)) = 48.888 fs), add the line TIMEUNIT INTERNAL
 - for more MD steps, change the MDSTEPS value
 - for a different temperature, change the line TEMPERATURE 0.616 (the units are kcal/mol!)

