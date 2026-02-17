MODULE UTILS_IO
   IMPLICIT NONE
   CONTAINS
      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)     
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
   
         FILE_UNIT = GETUNIT()
         IF (APPEND) THEN
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN', POSITION='APPEND')
         ELSE
            OPEN(UNIT=FILE_UNIT, FILE=FILE_NAME, STATUS='UNKNOWN')
         END IF   
      END SUBROUTINE FILE_OPEN

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

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
         INTEGER :: UNITNUM
         ! start checking for available units > 103, to avoid system default units
         ! 100, 101 and 102 are stdin, stdout and stderr respectively.
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
      
      !> Routine to parse a line into a list of words
      !>
      !> @param[in] LINE - input line to be aprsed
      !> @param[in] NWORDS - Number of words (for size of output array)
      !> @param[out] WORDSOUT - List of words from input line
      SUBROUTINE READLINE(LINE,NWORDS,WORDSOUT)
         CHARACTER(*), INTENT(IN) :: LINE
         INTEGER, INTENT(IN) :: NWORDS
         CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
         INTEGER:: J1,START_IND,END_IND,J2
         CHARACTER(35) :: WORD
         START_IND=0
         END_IND=0
         J1=1
         J2=0
         DO WHILE(J1.LE.LEN(LINE))
            IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
               START_IND=J1
            ENDIF
            IF (START_IND.GT.0) THEN
               IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
               IF (J1.EQ.LEN(LINE)) END_IND=J1
               IF (END_IND.GT.0) THEN
                  J2=J2+1
                  WORD=LINE(START_IND:END_IND)
                  WORDSOUT(J2)=TRIM(WORD)
                  START_IND=0
                  END_IND=0
               ENDIF
            ENDIF
            J1=J1+1
         ENDDO
      END SUBROUTINE READLINE
END MODULE UTILS_IO