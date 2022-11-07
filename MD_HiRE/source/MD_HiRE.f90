! MD engine for HiRE

!> @author 
!> Dr Konstantin Roeder, University of Cambridge, 2022
!>
!> @file
!> Main program to run MD simulation with HiRE force field
PROGRAM MD_HIRE
   USE NUMKIND
   USE MD_COMMONS, ONLY: MYUNIT
   USE MD_SETUP, ONLY: READ_SETTINGS, SETUP_POTENTIAL, START_TRACKING
   USE MD_SIMULATION
   USE MD_UTILS, ONLY: REPORT_PARAMS, MD_START
   IMPLICIT NONE
   REAL(KIND=REAL64) :: TSTART, TEND, TSETUP, TRUNS

   CALL CPU_TIME(TSTART)
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
   WRITE(MYUNIT,'(A)') " mdhire> Calling velocity initialisation"
   CALL INITIALISE_VEL()
   CALL ZERO_STEP()
   CALL CPU_TIME(TSETUP)

   ! 3. Run MD steps
   CALL RUN_MD()
   CALL CPU_TIME(TRUNS)

   ! 4. Finish run
   CALL MD_FINISH()

   CALL CPU_TIME(TEND)

   WRITE(MYUNIT,'(A)') " "  
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Total time:      ", TEND-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Setup time:      ", TSETUP-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Simulation time: ", TRUNS-TSETUP   
   CLOSE(MYUNIT) 
END PROGRAM MD_HIRE
