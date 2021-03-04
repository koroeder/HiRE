     
      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
         INTEGER                       :: GETUNIT
         
         FILE_UNIT = GETUNIT()
         IF (APPEND) THEN
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN', POSITION='APPEND')
         ELSE
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN')
         END IF
      
      END SUBROUTINE FILE_OPEN

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
      !
      ! start checking for available units > 103, to avoid system default units
      ! 100, 101 and 102 are stdin, stdout and stderr respectively.
      ! 
         INTEGER :: UNITNUM

         INUSE=.TRUE.
         UNITNUM=103

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT

      INTEGER FUNCTION FILE_LENGTH(FILE_NAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         INTEGER                       :: FILE_UNIT
         INTEGER                       :: IO_STATUS

         CALL FILE_OPEN(FILE_NAME, FILE_UNIT, .FALSE.)
         FILE_LENGTH = 0
         DO
            READ(FILE_UNIT, *, IOSTAT=IO_STATUS)
            IF (IO_STATUS /= 0) EXIT
            FILE_LENGTH = FILE_LENGTH + 1
         ENDDO

         CLOSE(FILE_UNIT)
      END FUNCTION FILE_LENGTH 

      SUBROUTINE VECNORM(VECTOR,NOPT)
      !
      ! normalises VECTOR (adapted from OPTIM)
      !
         USE PREC
         IMPLICIT NONE
         INTEGER :: J1 
         REAL(KIND = REAL64) :: DUMMY 
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND = REAL64), INTENT(INOUT) :: VECTOR(NOPT)

         DUMMY=0.0D0
         DO J1=1,NOPT
            DUMMY=DUMMY+VECTOR(J1)**2
         ENDDO

         IF (DUMMY.GT.0.0D0) THEN
            DUMMY=1.0D0/DSQRT(DUMMY)
            DO J1=1,NOPT
               VECTOR(J1)=VECTOR(J1)*DUMMY
            ENDDO
         ELSE
            PRINT *,'WARNING: zero size vector passed to VECNORM - STOP'
            STOP
         ENDIF

      END SUBROUTINE VECNORM

