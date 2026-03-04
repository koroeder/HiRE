PROGRAM HIRE_CALC
   USE HIRE_INTERFACE
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE
   
   INTEGER :: NARGS, NATOMS
   INTEGER :: I
!   INTEGER, PARAMETER :: REAL64 = SELECTED_REAL_KIND(15, 307)
   INTEGER, PARAMETER :: XUNIT = 11 
   INTEGER, PARAMETER :: STDOUT = 6 
   CHARACTER(LEN=50) :: TOPNAME !Name of topology file 
   CHARACTER(LEN=50) :: SCALEDATNAME !Name of scale.dat file 
   REAL(KIND = REAL64), ALLOCATABLE :: X(:), GRAD(:)
   REAL(KIND = REAL64) :: E, RMS
   
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
   ALLOCATE(X(3*NATOMS),GRAD(3*NATOMS))
   OPEN(XUNIT, FILE="start")
   READ(XUNIT, *) (X(I), I=1,3*NATOMS)
   CLOSE(XUNIT)
   WRITE(STDOUT,'(A)') "Read coordinates"
   ! Calculate energy and gradient
#ifdef FOR_PROFILING 
   DO I=1,500
#endif
      CALL HIRE_ENERGY_GRAD(3*NATOMS, X, E, GRAD, .FALSE.)
#ifdef FOR_PROFILING
   END DO
#endif
   WRITE(STDOUT,'(A,F15.6)') "Energy:   ", E
   RMS=MAX(SQRT(SUM(GRAD(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100 )
   WRITE(STDOUT,'(A,G15.6)') "RMS force: ", RMS
   CALL HIRE_DECOMPE(STDOUT)
   FLUSH(STDOUT)
   ! save a pdb file with chainid
   CALL DUMP_PDB(3*NATOMS, X, "conf_chainid.pdb",.TRUE.)
   ! Finish up
   CALL TERMINATE_HIRE()
END PROGRAM HIRE_CALC
