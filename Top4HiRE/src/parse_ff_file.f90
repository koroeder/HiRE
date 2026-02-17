MODULE PARSE_FF
   USE PREC_HIRE
   USE FF_GLOBALS

   CONTAINS
      SUBROUTINE PARSE_FF_FILES(NFF,FFFILES)
         USE UTILS_IO, ONLY: FILE_OPEN
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NFF
         CHARACTER(LEN=30), INTENT(IN) :: FFFILES(NFF)
         INTEGER :: IOSTAT, INUNIT
         CHARACTER(LEN=200) :: LINE
         INTEGER :: J,N1,N2,N3,N4,N5,N6,N7,N8,N9

         ! allocate variables by checking all FF files
         DO J=1,NFF
            CALL FILE_OPEN(FFFILES(J), INUNIT, .FALSE.)
            IOSTAT=0
            DO WHILE (IOSTAT.EQ.0)
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
               IF (TRIM(ADJUSTL(LINE)).EQ."#FILECONTENT") THEN
                  READ(LINE,*) N1,N2,N3,N4,N5,N6,N7,N8,N9
                  EXIT
               END IF
            END DO
            CLOSE(INUNIT)
            NBINTER = NBINTER + N1
            NBINTRA = NBINTRA + N2
            NAINTER = NAINTER + N3
            NAINTRA = NAINTRA + N4
            NQINTER = NQINTER + N5
            NQINTRA = NQINTRA + N6
            NDINTER = NDINTER + N7
            NDINTRA = NDINTRA + N8            
         END DO
         CALL ALLOCATE_FF_GLOBALS()

         N1=0;N2=0;N3=0;N4=0;N5=0;N6=0;N7=0;N8=0
         DO J=1,NFF
            CALL READ_FF_FILE(FFFILES(J),N1,N2,N3,N4,N5,N6,N7,N8)
         END DO
      END SUBROUTINE PARSE_FF_FILES

      SUBROUTINE READ_FF_FILE(FNAME,N1CURR,N2CURR,N3CURR,N4CURR,N5CURR,N6CURR,N7CURR,N8CURR)
         USE UTILS_IO, ONLY: FILE_OPEN
         IMPLICIT NONE
         CHARACTER(LEN=30), INTENT(IN) :: FNAME
         INTEGER, INTENT(INOUT) :: N1CURR,N2CURR,N3CURR,N4CURR,N5CURR,N6CURR,N7CURR,N8CURR
         INTEGER :: INUNIT
         CHARACTER(LEN=200) :: LINE
         INTEGER :: IOSTAT
         INTEGER :: N1,N2,N3,N4,N5,N6,N7,N8,N9

         CALL FILE_OPEN(FNAME, INUNIT, .FALSE.)

         IOSTAT=0
         DO WHILE (IOSTAT.EQ.0)
            READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
            IF (TRIM(ADJUSTL(LINE)).EQ."#FILECONTENT") THEN
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
               READ(LINE,*) N1,N2,N3,N4,N5,N6,N7,N8,N9
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#BONDS")
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N1+N2

               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#BONDS")
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N3+N4

               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#BONDS")
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N5+N6

               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#BONDS")
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N9

               END DO            
            END IF
         END DO
         CLOSE(INUNIT)
      END SUBROUTINE READ_FF_FILE


END MODULE PARSE_FF