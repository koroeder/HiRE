MODULE MD_UTILS
   CONTAINS
      SUBROUTINE ALLOC_COMMONS()
         USE MD_COMMONS, ONLY: NATOMS, ACC, VEL, MASSES, COORDS
         IMPLICIT NONE
         CALL DEALLOC_COMMONS()
         ALLOCATE(COORDS(3*NATOMS))
         ALLOCATE(ACC(3*NATOMS))
         ALLOCATE(VEL(3*NATOMS))
         ALLOCATE(MASSES(NATOMS))
      END SUBROUTINE ALLOC_COMMONS

      SUBROUTINE DEALLOC_COMMONS()
         USE MD_COMMONS, ONLY: ACC, VEL, COORDS, MASSES
         IMPLICIT NONE        
         IF (ALLOCATED(COORDS)) DEALLOCATE(COORDS)
         IF (ALLOCATED(ACC)) DEALLOCATE(ACC)
         IF (ALLOCATED(VEL)) DEALLOCATE(VEL)
         IF (ALLOCATED(MASSES)) DEALLOCATE(MASSES)
      END SUBROUTINE DEALLOC_COMMONS    
      
      SUBROUTINE TERMINATE_ERR(HIREINIT, ALLOCT)
         USE HIRE_INTERFACE, ONLY: TERMINATE_HIRE
         USE MD_COMMONS, ONLY: MYUNIT
         IMPLICIT NONE
         LOGICAL, INTENT(IN) :: HIREINIT, ALLOCT
         WRITE(MYUNIT,'(A)') " terminate_err> Terminate simulation due to error"
         IF (HIREINIT) CALL TERMINATE_HIRE()
         IF (ALLOCT) CALL DEALLOC_COMMONS()
         CLOSE(MYUNIT)
         STOP
      END SUBROUTINE TERMINATE_ERR

      SUBROUTINE REPORT_PARAMS()
         USE MD_COMMONS
         IMPLICIT NONE
         WRITE(MYUNIT,'(A)') " Molecular dynamics simulation for HiRE "
         WRITE(MYUNIT,'(A)') " ______________________________________ "
         WRITE(MYUNIT,'(A)') " "
         WRITE(MYUNIT,'(A,I10,A)') " settings> Run MD simulation for ", MDSTEPS, " steps"
         WRITE(MYUNIT,'(A,F6.2)') " settings> Time step for simulation:          ", DT
         IF (MDMETHOD.EQ."VV") THEN
            WRITE(MYUNIT,'(A)') " settings> MD simulation will use Velocity-Verlet"
         ELSE IF (MDMETHOD.EQ."LD") THEN
            WRITE(MYUNIT,'(A)') " settings> MD simulation will use Langevin dynamics"
            WRITE(MYUNIT,'(A,F6.2)') " settings> Gamma value for Langevin dynamics: ", GAMMA
         ELSE 
            WRITE(MYUNIT,'(2A)') " settings> MD method not recognised: ", MDMETHOD
            CALL TERMINATE_ERR(.FALSE., .FALSE.)
         END IF
         WRITE(MYUNIT,'(A,F8.2)') " settings> Temperature for MD simulation:     ", TEMP
         WRITE(MYUNIT,'(A)') " "
      END SUBROUTINE REPORT_PARAMS

      SUBROUTINE MD_START()
         USE MD_COMMONS, ONLY: MYUNIT, NTASKS, TASKID
         USE FILE_UTILS, ONLY: FILE_OPEN, FILE_EXIST
         IMPLICIT NONE
         CHARACTER(LEN=6) :: STRID

         IF (FILE_EXIST("mddata")) THEN
            IF (NTASKS.EQ.1) THEN
               CALL FILE_OPEN("mdout.log",MYUNIT,.TRUE.)
            ELSE
               
               WRITE(STRID,'(I6)') TASKID
               CALL FILE_OPEN("mdout.log."//ADJUSTL(TRIM(STRID)),MYUNIT,.TRUE.)
            END IF
         ELSE
            WRITE(*,'(A)') " Cannot locate input file"
            STOP
         END IF
      END SUBROUTINE MD_START

      SUBROUTINE RUNMIN(X)
         USE NUMKIND
         USE MD_COMMONS, ONLY: NATOMS, MYUNIT
         USE MINIMISATION, ONLY: MINIMISE, COLDFUSION
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: X(3*NATOMS)
         REAL(KIND=REAL64) :: ENERGY
         INTEGER :: ITDONE
         LOGICAL :: MFLAG
         
         ENERGY=0.0D0

         CALL MINIMISE(3*NATOMS,X,ENERGY,ITDONE,MFLAG,MYUNIT)
         IF (COLDFUSION) THEN
            WRITE(MYUNIT,'(A)') " runmin> Cold fusion occured, minimisatio failed - STOP"
            CALL TERMINATE_ERR(.TRUE., .TRUE.)
         ELSE IF (.NOT.MFLAG) THEN
            WRITE(MYUNIT,'(A)') " runmin> Minimisation did not converge - STOP"
            CALL TERMINATE_ERR(.TRUE., .TRUE.)
         ELSE
            WRITE(MYUNIT,*) " runmin> Minimisation converged in ", ITDONE, " steps"
            WRITE(MYUNIT,*) "         Energy of minimum: ", ENERGY
         END IF
         WRITE(MYUNIT,*) ""
      END SUBROUTINE RUNMIN

      SUBROUTINE SET_DERIVED_PARAMS()
         USE MD_COMMONS, ONLY: DT, HDT, GAMMA, GFRIC, NATOMS, NOPT
         IMPLICIT NONE
         HDT = 0.5*DT
         GFRIC = 1.0D0 - GAMMA*HDT
      END SUBROUTINE SET_DERIVED_PARAMS
END MODULE MD_UTILS
