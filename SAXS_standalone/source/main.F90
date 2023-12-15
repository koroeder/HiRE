PROGRAM HIRE_CALC
   USE HIRE_INTERFACE
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE
   
   INTEGER :: NARGS, NATOMS
   INTEGER :: I
!   INTEGER, PARAMETER :: REAL64 = SELECTED_REAL_KIND(15, 307)
   INTEGER, PARAMETER :: XUNIT = 11 
   INTEGER, PARAMETER :: STDOUT = 6 
   INTEGER, PARAMETER :: PUNIT = 12
   CHARACTER(LEN=50) :: TOPNAME !Name of topology file 
   CHARACTER(LEN=50) :: SCALEDATNAME !Name of scale.dat file 
   REAL(KIND = REAL64), ALLOCATABLE :: X(:), GRAD(:), SAXSFORCE(:)
   REAL(KIND = REAL64) :: E, RMS, ESAXS
   CHARACTER(LEN=12) :: DUMMY
   LOGICAL :: YESNO
   LOGICAL :: CALC_HYDRATIONT = .FALSE. !calculate hydration
   LOGICAL :: CALC_FORCET = .TRUE. !calculate SAXS force
   LOGICAL           :: SAXSMODULT = .FALSE.      !modulate SAXS force with decreasing and periodic factor?
   REAL(KIND=REAL64) :: SAXSINVSIG = 2.5          !broadness of the periodic modulation
   REAL(KIND=REAL64) :: SAXSWAVE = 100.0          !length of the decreasing modulation 
   REAL(KIND=REAL64) :: SAXSOFFI = 0.0D0          !step variable for periodic modulation
   REAL(KIND=REAL64) :: SAXSMODI = 1.0D0          !step variable for decreasing modulation
   LOGICAL           :: SAXSPRINT = .FALSE.       !print curves and scores to unit
   LOGICAL           :: SAXSSAVET  = .FALSE.      !save computed SAXS at MC transition?
   REAL(KIND=REAL64) :: SAXSMAX = 0.35             !max value for SAXS vector
   LOGICAL           :: SAXSSOLT = .FALSE.        !vacuum (false) or solution (true) calculation
   LOGICAL           :: REFINET = .FALSE.         !refine hydration
   INTEGER           :: NWATLAY = 1               !number of hydration layers
   REAL(KIND=REAL64) :: WATRAD = 2.55             !water radius
   REAL(KIND=REAL64) :: WATW = 0.037              !water contrast
   ! check number of arguments
   NARGS = COMMAND_ARGUMENT_COUNT()
   ! We expect two arguments, the topology file and the scale.dat file
   IF (NARGS.EQ.2) THEN
      CALL GET_COMMAND_ARGUMENT(1, TOPNAME)
      CALL GET_COMMAND_ARGUMENT(2, SCALEDATNAME) 
      WRITE(STDOUT,'(A)') "Call setup for HiRE"
      CALL HIRE_INITIALISE(TOPNAME, SCALEDATNAME, NATOMS)
      WRITE(STDOUT,'(A)') "Finished setup"
   ELSE
      WRITE(STDOUT,'(A,I4)') "Expecting two arguments, but got ", NARGS
      STOP
   END IF
   ! Read in coordinates
   ALLOCATE(X(3*NATOMS),GRAD(3*NATOMS),SAXSFORCE(3*NATOMS))
   OPEN(XUNIT, FILE="start")
   READ(XUNIT, *) (X(I), I=1,3*NATOMS)
   CLOSE(XUNIT)
   WRITE(STDOUT,'(A)') "Read coordinates"
   ! Calculate energy and gradient as reference
   CALL HIRE_ENERGY_GRAD(3*NATOMS, X, E, GRAD, .FALSE.)
   WRITE(STDOUT,'(A,F15.6)') "Energy without SAXS:   ", E
   RMS=MAX(SQRT(SUM(GRAD(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100 )
   WRITE(STDOUT,'(A,G15.6)') "RMS force without SAXS: ", RMS

   !read in SAXS parameters
   INQUIRE(FILE="parameters.SAXS", EXIST=YESNO)
   IF (YESNO) THEN
      WRITE(STDOUT, '(A)') " Reading SAXS parameters from file parameters.SAXS"
      OPEN(PUNIT, FILE="parameters.SAXS")
      READ(PUNIT, *) DUMMY, SAXSMODULT
      READ(PUNIT, *) DUMMY, SAXSINVSIG
      READ(PUNIT, *) DUMMY, SAXSMAX
      READ(PUNIT, *) DUMMY, SAXSSOLT
      READ(PUNIT, *) DUMMY, REFINET
      READ(PUNIT, *) DUMMY, WATRAD
      READ(PUNIT, *) DUMMY, NWATLAY
      READ(PUNIT, *) DUMMY, CALC_HYDRATIONT
      READ(PUNIT, *) DUMMY, CALC_FORCET
      CLOSE(PUNIT)
   ELSE
      WRITE(STDOUT, '(A)') " Using default SAXS parameters"
   END IF
   FLUSH(STDOUT)
   !setup SAXS
   CALL SETUP_SAXS(.TRUE.,SAXSPRINT,SAXSMODULT,SAXSINVSIG,SAXSMAX, &
   SAXSSOLT,REFINET,WATRAD,NWATLAY)
   !calcualte SAXS energy and force
#ifdef FOR_PROFILING 
   DO I=1,10
#endif
      CALL HIRE_SAXS_FORCE(3*NATOMS,X,ESAXS,SAXSFORCE,CALC_FORCET,CALC_HYDRATIONT)
#ifdef FOR_PROFILING
   END DO
#endif
   !CALL HIRE_SAXS_FORCE(3*NATOMS,X,ESAXS,SAXSFORCE,1)
   !WRITE(STDOUT,'(A,F15.6)') "SAXS Energy:   ", ESAXS
   !RMS=MAX(SQRT(SUM(SAXSFORCE(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100 )
   !WRITE(STDOUT,'(A,G15.6)') "SAXS RMS force: ", RMS
   !FLUSH(STDOUT)   
   WRITE(STDOUT,'(A,F15.6)') "SAXS Energy:   ", ESAXS
   RMS=MAX(SQRT(SUM(SAXSFORCE(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100 )
   WRITE(STDOUT,'(A,G15.6)') "SAXS RMS force: ", RMS
   FLUSH(STDOUT)  
   ! Finish up
   CALL TERMINATE_HIRE()
END PROGRAM HIRE_CALC
