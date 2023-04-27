! MD engine for HiRE

!> @author 
!> Dr Konstantin Roeder, University of Cambridge, 2022
!>
!> @file
!> Main program to run MD simulation with HiRE force field
PROGRAM MD_HIRE
   USE NUMKIND
   USE MD_COMMONS, ONLY: MYUNIT, THERMINIT, TASKID, NTASKS
   USE MD_SETUP, ONLY: READ_SETTINGS, SETUP_POTENTIAL, START_TRACKING
   USE MD_SIMULATION
   USE MD_UTILS, ONLY: REPORT_PARAMS, MD_START, SEED_RANDOM, TEST_RANDOM
#ifdef MPI
   USE MPI_UTILS, ONLY: COMMUNICATE_SETTINGS, REPORT_PARAMS_MPI, START_TRACKING_MPI
#endif
   IMPLICIT NONE
   REAL(KIND=REAL64) :: TSTART, TEND, TSETUP, TEQ, TRUNS

#ifdef MPI
   INCLUDE 'mpif.h'
   INTEGER :: ERROR

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
   ! b) Read in all the simulation settings
   WRITE(MYUNIT,'(A)') " mdhire> Starting setup"
   CALL READ_SETTINGS()
   CALL SEED_RANDOM()
   ! communicate the variables
#ifdef MPI
   CALL COMMUNICATE_SETTINGS()
#endif
   ! c) Report the settings for the simulation to the output file
#ifdef MPI
   CALL REPORT_PARAMS_MPI()
#else
   CALL REPORT_PARAMS()
#endif
   ! d) Initialise potential and get coordinates, mass and names from HiRE
   CALL SETUP_POTENTIAL()
   ! e) Open tracking files
#ifdef MPI
   CALL START_TRACKING_MPI
#else
   CALL START_TRACKING()
#endif
   CALL CPU_TIME(TSETUP)
   WRITE(MYUNIT,'(A)') " mdhire> Completed setup - starting initialisation"
   ! 2. Initialise velocities for simulation and get initial state
   CALL ZERO_STEP()
   CALL CPU_TIME(TEQ)
#ifdef MPI
   CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR)
#endif
   ! 3. Run MD steps
   WRITE(MYUNIT,'(A)') " mdhire> Starting MD simulation"
   CALL RUN_MD()
   CALL CPU_TIME(TRUNS)

#ifdef MPI
   CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR)
#endif
   WRITE(MYUNIT,'(A)') " mdhire> Completed simulation - terminating now"   
   ! 4. Finish run
   CALL MD_FINISH()

   CALL CPU_TIME(TEND)

   WRITE(MYUNIT,'(A)') " "  
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Total time:      ", TEND-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Setup time:      ", TSETUP-TSTART
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Thermalisation : ", TEQ-TSETUP 
   WRITE(MYUNIT,'(A,F10.2)') " mdhire> Simulation time: ", TRUNS-TEQ  
   CLOSE(MYUNIT) 
#ifdef MPI
   CALL MPI_BARRIER(MPI_COMM_WORLD,ERROR)
   CALL MPI_FINALIZE(ERROR)
#endif
END PROGRAM MD_HIRE
