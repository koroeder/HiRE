MODULE PARSE_FF
   USE PREC_HIRE
   USE FF_GLOBALS

   CONTAINS
      SUBROUTINE PARSE_FF_FILES(NFF,FFFILES)
         USE UTILS_IO, ONLY: FILE_OPEN
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NFF
         CHARACTER(LEN=256), INTENT(IN) :: FFFILES(NFF)
         INTEGER :: IOSTAT, INUNIT
         CHARACTER(LEN=200) :: LINE
         INTEGER :: I,J,N1,N2,N3,N4,N5,N6,N7,N8,N9

         ! allocate variables by checking all FF files
         DO J=1,NFF
            CALL FILE_OPEN(FFFILES(J), INUNIT, .FALSE.)
            IOSTAT=0
            DO WHILE (IOSTAT.EQ.0)
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
               IF (TRIM(ADJUSTL(LINE)).EQ."#FILECONTENT") THEN
                  READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
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
            WRITE(*,*) " Using force field data file: ", FFFILES(J)
         END DO
      END SUBROUTINE PARSE_FF_FILES

      SUBROUTINE READ_FF_FILE(FNAME,N1CURR,N2CURR,N3CURR,N4CURR,N5CURR,N6CURR,N7CURR,N8CURR)
         USE UTILS_IO, ONLY: FILE_OPEN, READLINE
         IMPLICIT NONE
         CHARACTER(LEN=256), INTENT(IN) :: FNAME
         INTEGER, INTENT(INOUT) :: N1CURR,N2CURR,N3CURR,N4CURR,N5CURR,N6CURR,N7CURR,N8CURR
         INTEGER :: INUNIT
         CHARACTER(LEN=200) :: LINE
         INTEGER, PARAMETER :: NENTRIES = 15
         CHARACTER(LEN=20) :: ENTRIES(NENTRIES)
         INTEGER :: IOSTAT
         INTEGER :: N1,N2,N3,N4,N5,N6,N7,N8,N9
         INTEGER :: INTDUMMY, INTTYPE, I, DIDX, TERMIDX
         INTEGER, ALLOCATABLE :: CURRTERMS(:), DIHTYPE(:)

         CALL FILE_OPEN(FNAME, INUNIT, .FALSE.)

         IOSTAT=0
         DO WHILE (IOSTAT.EQ.0)
            READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
            IF (TRIM(ADJUSTL(LINE)).EQ."#FILECONTENT") THEN
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
               READ(LINE,*) N1,N2,N3,N4,N5,N6,N7,N8,N9
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#BONDS") THEN
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N1+N2
                  READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
                  CALL READLINE(LINE,NENTRIES,ENTRIES)
                  READ(ENTRIES(6),'(I4)') INTDUMMY
                  IF (INTDUMMY.GE.1) THEN !internulceotide
                     N1CURR = N1CURR + 1
                     BINTER(N1CURR)%AT1 = ENTRIES(2)
                     BINTER(N1CURR)%AT2 = ENTRIES(3)
                     READ(ENTRIES(4),'(F12.4)') BINTER(N1CURR)%KSPR
                     READ(ENTRIES(5),'(F12.4)') BINTER(N1CURR)%REQ
                     READ(ENTRIES(1),'(I4)') BINTERTYPE(N1CURR)
                  ELSE
                     N2CURR = N2CURR + 1
                     BINTRA(N2CURR)%AT1 = ENTRIES(2)
                     BINTRA(N2CURR)%AT2 = ENTRIES(3)
                     READ(ENTRIES(4),'(F12.4)') BINTRA(N2CURR)%KSPR
                     READ(ENTRIES(5),'(F12.4)') BINTRA(N2CURR)%REQ
                     READ(ENTRIES(1),'(I4)') BINTRATYPE(N2CURR)
                  END IF
               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#ANGLES") THEN
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N3+N4
                  READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
                  CALL READLINE(LINE,NENTRIES,ENTRIES)
                  READ(ENTRIES(7),'(I4)') INTDUMMY
                  IF (INTDUMMY.GE.1) THEN !internulceotide
                     N3CURR = N3CURR + 1
                     AINTER(N3CURR)%AT1 = ENTRIES(2)
                     AINTER(N3CURR)%AT2 = ENTRIES(3)
                     AINTER(N3CURR)%AT3 = ENTRIES(4)
                     READ(ENTRIES(5),'(F12.4)') AINTER(N3CURR)%KSPR
                     READ(ENTRIES(6),'(F12.4)') AINTER(N3CURR)%TEQ
                     READ(ENTRIES(1),'(I4)') AINTERTYPE(N3CURR)
                  ELSE
                     N4CURR = N4CURR + 1
                     AINTRA(N4CURR)%AT1 = ENTRIES(2)
                     AINTRA(N4CURR)%AT2 = ENTRIES(3)
                     AINTRA(N4CURR)%AT3 = ENTRIES(4)
                     READ(ENTRIES(5),'(F12.4)') AINTRA(N4CURR)%KSPR
                     READ(ENTRIES(6),'(F12.4)') AINTRA(N4CURR)%TEQ
                     READ(ENTRIES(1),'(I4)') AINTRATYPE(N4CURR)
                  END IF
               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#QANGLES") THEN
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N5+N6
                  READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
                  CALL READLINE(LINE,NENTRIES,ENTRIES)
                  READ(ENTRIES(10),'(I4)') INTDUMMY
                  IF (INTDUMMY.GE.1) THEN !internulceotide
                     N5CURR = N5CURR + 1
                     QINTER(N5CURR)%AT1 = ENTRIES(2)
                     QINTER(N5CURR)%AT2 = ENTRIES(3)
                     QINTER(N5CURR)%AT3 = ENTRIES(4)
                     READ(ENTRIES(5),'(F12.4)') QINTER(N5CURR)%TTA
                     READ(ENTRIES(6),'(F12.4)') QINTER(N5CURR)%A1
                     READ(ENTRIES(7),'(F12.4)') QINTER(N5CURR)%A2
                     READ(ENTRIES(8),'(F12.4)') QINTER(N5CURR)%A3
                     READ(ENTRIES(9),'(F12.4)') QINTER(N5CURR)%A5
                     READ(ENTRIES(1),'(I4)') QINTERTYPE(N5CURR)
                  ELSE
                     N6CURR = N6CURR + 1
                     QINTRA(N6CURR)%AT1 = ENTRIES(2)
                     QINTRA(N6CURR)%AT2 = ENTRIES(3)
                     QINTRA(N6CURR)%AT3 = ENTRIES(4)
                     READ(ENTRIES(5),'(F12.4)') QINTRA(N6CURR)%TTA
                     READ(ENTRIES(6),'(F12.4)') QINTRA(N6CURR)%A1
                     READ(ENTRIES(7),'(F12.4)') QINTRA(N6CURR)%A2
                     READ(ENTRIES(8),'(F12.4)') QINTRA(N6CURR)%A3
                     READ(ENTRIES(9),'(F12.4)') QINTRA(N6CURR)%A5
                     READ(ENTRIES(1),'(I4)') QINTRATYPE(N6CURR)
                  END IF
               END DO
            ELSE IF (TRIM(ADJUSTL(LINE)).EQ."#DIHEDRALS") THEN
               !allocate helper variables
               ALLOCATE(CURRTERMS(N7+N8),DIHTYPE(N7+N8))
               CURRTERMS(1:N7+N8) = 0
               DIHTYPE(1:N7+N8) = -1
               READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE !comment line
               DO I=1,N9
                  READ(INUNIT,'(A)',IOSTAT=IOSTAT) LINE
                  CALL READLINE(LINE,NENTRIES,ENTRIES)
                  READ(ENTRIES(11),'(I4)') INTDUMMY
                  READ(ENTRIES(2),'(I4)') INTTYPE
                  IF (INTDUMMY.GE.1) THEN !internulceotide
                     !first determine whether this is a new dihedral, or whether we already have some terms
                     IF (DIHTYPE(INTTYPE).EQ.-1) THEN
                        N7CURR = N7CURR + 1
                        DIHTYPE(INTTYPE) = N7CURR
                        DIDX = DIHTYPE(INTTYPE)
                     ELSE
                        DIDX = DIHTYPE(INTTYPE)
                     END IF
                     CURRTERMS(INTTYPE) = CURRTERMS(INTTYPE) + 1
                     TERMIDX = CURRTERMS(INTTYPE)
                     IF (TERMIDX.GT.MAXNTERM) THEN
                        WRITE(*,*) "Error: too many terms for dihedral type ", INTTYPE
                        STOP
                     END IF
                     DINTER(DIDX,TERMIDX)%AT1 = ENTRIES(3)
                     DINTER(DIDX,TERMIDX)%AT2 = ENTRIES(4)
                     DINTER(DIDX,TERMIDX)%AT3 = ENTRIES(5)
                     DINTER(DIDX,TERMIDX)%AT4 = ENTRIES(6)
                     READ(ENTRIES(7),'(F12.4)') DINTER(DIDX,TERMIDX)%K
                     READ(ENTRIES(8),'(F12.4)') DINTER(DIDX,TERMIDX)%PHASE
                     READ(ENTRIES(9),'(F8.0)') DINTER(DIDX,TERMIDX)%PERIOD
                     READ(ENTRIES(10),'(F12.4)') DINTER(DIDX,TERMIDX)%YOFF
                     READ(ENTRIES(1),'(I4)') DINTERTYPE(DIDX)
                  ELSE
                     !first determine whether this is a new dihedral, or whether we already have some terms
                     IF (DIHTYPE(INTTYPE).EQ.-1) THEN
                        N8CURR = N8CURR + 1
                        DIHTYPE(INTTYPE) = N8CURR
                        DIDX = DIHTYPE(INTTYPE)
                     ELSE
                        DIDX = DIHTYPE(INTTYPE)
                     END IF
                     CURRTERMS(INTTYPE) = CURRTERMS(INTTYPE) + 1
                     TERMIDX = CURRTERMS(INTTYPE)
                     IF (TERMIDX.GT.MAXNTERM) THEN
                        WRITE(*,*) "Error: too many terms for dihedral type ", INTTYPE
                        STOP
                     END IF
                     DINTRA(DIDX,TERMIDX)%AT1 = ENTRIES(3)
                     DINTRA(DIDX,TERMIDX)%AT2 = ENTRIES(4)
                     DINTRA(DIDX,TERMIDX)%AT3 = ENTRIES(5)
                     DINTRA(DIDX,TERMIDX)%AT4 = ENTRIES(6)
                     READ(ENTRIES(7),'(F12.4)') DINTRA(DIDX,TERMIDX)%K
                     READ(ENTRIES(8),'(F12.4)') DINTRA(DIDX,TERMIDX)%PHASE
                     READ(ENTRIES(9),'(F8.0)') DINTRA(DIDX,TERMIDX)%PERIOD
                     READ(ENTRIES(10),'(F12.4)') DINTRA(DIDX,TERMIDX)%YOFF
                     READ(ENTRIES(1),'(I4)') DINTRATYPE(DIDX)
                  END IF
               END DO
               DEALLOCATE(CURRTERMS,DIHTYPE)
            END IF
         END DO
         CLOSE(INUNIT)
      END SUBROUTINE READ_FF_FILE
END MODULE PARSE_FF
