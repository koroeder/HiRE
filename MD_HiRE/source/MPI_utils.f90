MODULE MD_MPI
   USE NUMKIND
   IMPLICIT NONE
   
   INTEGER :: ERR_CODE_MPI
   INTEGER :: REMD_TAG = 12
   CONTAINS

      SUBROUTINE COMMUNICATE_SETTINGS()
         USE MD_COMMONS
         IMPLICIT NONE
         REAL(KIND = REAL64), ALLOCATABLE :: TARRAY(:), STEP
         INTEGER :: J
         
         
         IF(TASKID.EQ.0) THEN
            IF (NREPLICA.NE.NTASKS) THEN
               WRITE(MYUNIT,'(A)') " comm_settings> WARNING - number of replica and number of tasks differ, &
                                     resetting them to be equal by changing the number of replicas"
               NREPLICA = NTASKS
            END IF
         END IF
         ! use BCAST for globals

         ! use SEND/RECV for temperature and lambda
         ALLOCATE(TARRAY(NREPLICA))
         STEP = (HIGHR - LOWR)/DBLE(NREPLICA-1)
         DO J=0,NREPLICA-1
            TARRAY(J+1) = LOWR + J*STEP
         END DO
         IF (REXMODE.EQ."T") THEN
            IF (TASKID.EQ.0) THEN
               TEMP = TARRAY(1)
               DO J=2,NREPLICA
                  CALL MPI_SEND(TARRAY(J),1,MPI_DOUBLE,J-1,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
               END DO
            ELSE 
               CALL MPI_REC(TEMP,1,MPI_DOUBLE,0,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
            END IF
         END IF
         IF (REXMODE.EQ."H") THEN
            IF (TASKID.EQ.0) THEN
               LAMBDA = TARRAY(1)
               DO J=2,NREPLICA
                  CALL MPI_SEND(TARRAY(J),1,MPI_DOUBLE,J-1,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
               END DO
            ELSE 
               CALL MPI_REC(LAMBDA,1,MPI_DOUBLE,0,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
            END IF
            CALL MPI_BCAST(TEMP,1,MPI_DOUBLE,0,MPI_COMM_WORLD,ERR_CODE_MPI)
         END IF           
      END SUBROUTINE COMMUNICATE_SETTINGS


END MODULE MD_MPI