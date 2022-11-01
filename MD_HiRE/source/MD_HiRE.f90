! MD engine for HiRE

!> @author 
!> Dr Konstantin Roeder, University of Cambridge, 2022
!>
!> @file
!> Main program to run MD simulation with HiRE force field
PROGRAM MD_HIRE
   USE MD_SETUP, ONLY: READ_SETTINGS, SETUP_POTENTIAL
   USE MD_SIMULATION
   USE MD_FINAL
   USE MD_UTILS, ONLY: REPORT_PARAMS
   
   ! 1. Initialise everything appropriately
   CALL READ_SETTINGS()
   CALL SETUP_POTENTIAL()

   ! 2. Report the settings for the simulation to the output file

END PROGRAM MD_HIRE
