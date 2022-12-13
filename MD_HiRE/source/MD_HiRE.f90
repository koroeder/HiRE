! MD engine for HiRE

!> @author 
!> Dr Konstantin Roeder, University of Cambridge, 2022
!>
!> @file
!> Main program to run MD simulation with HiRE force field
PROGRAM MD_HIRE
   USE NUMKIND
   USE MD_COMMONS, ONLY: MYUNIT, THERMINIT
   USE MD_SETUP, ONLY: READ_SETTINGS, SETUP_POTENTIAL, START_TRACKING
   USE MD_SIMULATION
   USE MD_UTILS, ONLY: REPORT_PARAMS, MD_START
   USE MPI_UTILS, ONLY: COMMUNICATE_SETTINGS
   IMPLICIT NONE
   REAL(KIND=REAL64) :: TSTART, TEND, TSETUP, TEQ, TRUNS
   INTEGER :: ERROR
#ifdef MPI
   CALL MPI_INIT(ERROR)
   CALL MPI_COMM_SIZE(MPI_COMM_WORLD, NTASKS, ERROR)
   CALL MPI_COMM_RANK(MPI_COMM_WORLD, TASKID, ERROR)
#else
   NTASKS = 1
   TASKID = 0
#endif


   CALL CPU_TIME(TSTART)
   ! 1. Simulation setup
   ! a) check the parameter input exists
   CALL MD_START()
   ! b) Read in all the simulation settings if we are tasks 0
   IF (TASKID.EQ.0) THEN
      CALL READ_SETTINGS()
   END IF
   ! communicate the variables
#ifdef MPI
   CALL COMMUNICATE_SETTINGS()
#endif
   ! c) Report the settings for the simulation to the output file
   CALL REPORT_PARAMS()
   ! d) Initialise potential and get coordinates, mass and names from HiRE
   CALL SETUP_POTENTIAL()
   ! e) Open tracking files
   CALL START_TRACKING()
   CALL CPU_TIME(TSETUP)
   
   ! 2. Initialise velocities for simulation and get initial state
   CALL ZERO_STEP()
   CALL CPU_TIME(TEQ)

   ! 3. Run MD steps
   CALL RUN_MD()
   CALL CPU_TIME(TRUNS)

   ! 4. Finish run
   CALL MD_FINISH()

   CALL CPU_TIME(TEND)

   WRITE(MYUNIT,'(A)') " "  
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Total time:      ", TEND-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Setup time:      ", TSETUP-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Thermalisation : ", TEQ-TSETUP 
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Simulation time: ", TRUNS-TEQ  
   CLOSE(MYUNIT) 
END PROGRAM MD_HIRE
