! MD engine for HiRE

!> @author 
!> Dr Konstantin Roeder, University of Cambridge, 2022
!>
!> @file
!> Main program to run MD simulation with HiRE force field
PROGRAM MD_HIRE
   USE MD_COMMONS, ONLY: MYUNIT
   USE MD_SETUP, ONLY: READ_SETTINGS, SETUP_POTENTIAL
   USE MD_SIMULATION
   USE MD_FINAL
   USE MD_UTILS, ONLY: REPORT_PARAMS, MD_START, START_TRACKING
   
   ! 1. Simulation setup
   ! a) check the parameter input exists
   CALL MD_START()
   ! b) Read in all the simulation settings
   CALL READ_SETTINGS()
   ! c) Initialise potential and get coordinates, mass and names from HiRE
   CALL SETUP_POTENTIAL()
   ! c) Report the settings for the simulation to the output file
   CALL REPORT_PARAMS()
   ! d) Open tracking files
   CALL START_TRACKING()

   ! 2. Initialise velocities for simulation and get initial state
   WRITE(MYUNIT,'(A)') ""
   WRITE(MYUNIT,'(A)') " mdhire> Calling velocity initialisation"
   CALL INITIALISE_VEL()
   CALL ZERO_STEP()

   ! 3. Run MD steps
END PROGRAM MD_HIRE
