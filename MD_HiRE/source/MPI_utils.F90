MODULE MPI_UTILS
   USE NUMKIND
   IMPLICIT NONE

   INTEGER :: ERR_CODE_MPI
   INTEGER :: REMD_TAG = 12
   CONTAINS

      SUBROUTINE REPORT_PARAMS_MPI()
         USE MD_COMMONS
         USE MD_UTILS, ONLY: TERMINATE_ERR
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
         WRITE(MYUNIT,'(A,I6,A,I6)') " settings> This is a replica-exchange simulation with ", NREPLICA, &
                                     ". This is replica ", TASKID+1
         IF (REXMODE.EQ.'T') THEN
            WRITE(MYUNIT,'(A)') " settings> REX type: Temperature-REX"
         ELSE IF (REXMODE.EQ.'H') THEN
            WRITE(MYUNIT,'(A)') " settings> REX type: Hamiltonian-REX"
            WRITE(MYUNIT,'(A,F8.2)') " settings> Lambda is: ", LAMBDA
         END IF
         WRITE(MYUNIT,'(A,F8.2)') " settings> Temperature for MD simulation:     ", TEMP
         WRITE(MYUNIT,'(A)') " "
      END SUBROUTINE REPORT_PARAMS_MPI

      SUBROUTINE START_TRACKING_MPI()
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE MD_COMMONS
         USE EXCHANGES, ONLY: REXUNIT
         IMPLICIT NONE
         CHARACTER(LEN=6) :: IDSTR

         WRITE(IDSTR,'(I6)') TASKID + 1
         CALL FILE_OPEN("md_energy."//TRIM(ADJUSTL(IDSTR))//".log",EUNIT,.TRUE.)
         CALL FILE_OPEN("md_coords."//TRIM(ADJUSTL(IDSTR))//".xyz",XUNIT,.TRUE.)
         CALL FILE_OPEN("md_temp."//TRIM(ADJUSTL(IDSTR))//".log",TEMPUNIT,.TRUE.)
         IF (RMSDT) CALL FILE_OPEN("md_rmsd."//TRIM(ADJUSTL(IDSTR))//".log",RUNIT,.TRUE.)
         IF (TASKID.EQ.0) CALL FILE_OPEN("md_rexid.log",REXUNIT,.TRUE.)
      END SUBROUTINE START_TRACKING_MPI

      SUBROUTINE COMMUNICATE_SETTINGS()
         USE MD_COMMONS
         USE HIRE_INTERFACE, ONLY:  SET_UNIV_SCALING
         USE FILE_UTILS, ONLY: FILE_EXIST, FILE_OPEN
         USE MD_UTILS, ONLY: TERMINATE_ERR
         IMPLICIT NONE
         REAL(KIND = REAL64), ALLOCATABLE :: TARRAY(:), STEP
#ifdef MPI       
         INTEGER :: J
         INCLUDE 'mpif.h'
         INTEGER MPISTATUS(MPI_STATUS_SIZE)
         INTEGER :: ERR_CODE_MPI, REMD_TAG
         REAL(KIND=REAL64) :: KST
         INTEGER :: TEMPSUNIT

         IF(TASKID.EQ.0) THEN
            IF (NREPLICA.NE.NTASKS) THEN
               WRITE(MYUNIT,'(A)') " comm_settings> WARNING  number of replica and number of tasks differ, &
                                     resetting them to be equal by changing the number of replicas"
               NREPLICA = NTASKS
            END IF
            ALLOCATE(TARRAY(NREPLICA))
            
            IF (READTEMPS) THEN
               IF (.NOT.FILE_EXIST(TEMPSFILE)) THEN
                  WRITE(MYUNIT,'(A)') " comm_settings> Cannot locate file with REX temperatures - STOP"
                  CALL TERMINATE_ERR(.FALSE., .FALSE.)
               END IF
               CALL FILE_OPEN(TEMPSFILE,TEMPSUNIT,.FALSE.)
               DO J=1,NREPLICA
                  CALL INPUT(ENDT, TEMPSUNIT)
                  CALL READF(TARRAY(J))         
                  IF (ENDT) THEN
                     WRITE(MYUNIT,*) "End of file before all entries were read: ", TEMPSFILE, " - STOP"
                     CALL TERMINATE_ERR(.FALSE., .FALSE.)
                  END IF
               END DO
               CLOSE(TEMPSUNIT)
            ELSE 
               ! STEP = (HIGHR - LOWR)/DBLE(NREPLICA-1)
               KST = LOG(HIGHR/LOWR)/DBLE(NREPLICA-1)
               DO J=0,NREPLICA-1
                  ! TARRAY(J+1) = LOWR + J*STEP
                  TARRAY(J+1) = LOWR*EXP(KST*J)
               END DO
            END IF
         END IF

         ! use SEND/RECV for temperature and lambda
         IF (REXMODE.EQ."T") THEN
            IF (TASKID.EQ.0) THEN
               TEMP = TARRAY(1)
               DO J=2,NREPLICA
                  REMD_TAG = J
                  CALL MPI_SSEND(TARRAY(J),1,MPI_DOUBLE,J-1,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
               END DO
            ELSE
               REMD_TAG = TASKID + 1
               CALL MPI_RECV(TEMP,1,MPI_DOUBLE,0,REMD_TAG,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
            END IF
         END IF
         IF (REXMODE.EQ."H") THEN
            IF (TASKID.EQ.0) THEN
               LAMBDA = TARRAY(1)
               DO J=2,NREPLICA
                  REMD_TAG = J
                  CALL MPI_SEND(TARRAY(J),1,MPI_DOUBLE,J-1,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
               END DO
            ELSE
               REMD_TAG = TASKID + 1
               CALL MPI_RECV(LAMBDA,1,MPI_DOUBLE,0,REMD_TAG,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
            END IF
            CALL MPI_BCAST(TEMP,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ERR_CODE_MPI)
            CALL SET_UNIV_SCALING(LAMBDA)
         END IF
#endif
      END SUBROUTINE COMMUNICATE_SETTINGS

END MODULE MPI_UTILS
