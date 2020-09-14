MODULE INPUTmod

  USE PREC

  INTEGER :: NITEMS, NCURRENT
  CHARACTER (LEN=200) :: LINE
  CHARACTER(25), ALLOCATABLE :: WORDSOUT(:)  


  CONTAINS

    SUBROUTINE INPUT(ENDT,KEYWORD,IR,UPT)
       IMPLICIT NONE
       LOGICAL, INTENT(OUT)      :: ENDT     !end of file?
       CHARACTER(*), INTENT(OUT) :: KEYWORD  !first word in line
       INTEGER, INTENT(IN)       :: IR       !file unit
       LOGICAL, INTENT(IN)       :: UPT      !logical to capatilise KEYWORD 

       INTEGER :: NENTRIES
      
10     READ(IR,'(A)',END=20,ERR=25) LINE
       NITEMS = 0
       CALL READ_LINE(NENTRIES)
       IF (NENTRIES.EQ.0) THEN
!          WRITE(*,'(A)') 'input> Blank line in data file, reading next line'
          GOTO 10
       ELSE
          NITEMS = NENTRIES
          KEYWORD = WORDSOUT(1)
!          WRITE(*,'(3A,I8)') 'input> Number of options specified for keyword ', KEYWORD, ' is ', NITEMS
       ENDIF
       IF (UPT) CALL L2U(KEYWORD)
       ENDT=.FALSE.
       NCURRENT=1
       RETURN
20     ENDT=.TRUE.
       RETURN   
25     STOP 'input> Error while reading next line'
    END SUBROUTINE

    SUBROUTINE READA(A)
       IMPLICIT NONE
       CHARACTER(*), INTENT(OUT) :: A
       
       NCURRENT = NCURRENT + 1
       IF (NCURRENT.GT.NITEMS) THEN
          CALL REPORT('input> Too few arguments provided, use default!',.TRUE.)
          RETURN
       ENDIF
       A=WORDSOUT(NCURRENT)
       RETURN
    END SUBROUTINE READA

    SUBROUTINE READF(F)
       IMPLICIT NONE
       REAL(KIND=REAL64), INTENT(OUT) :: F
    
       NCURRENT = NCURRENT + 1
       IF (NCURRENT.GT.NITEMS) THEN
          CALL REPORT('input> Too few arguments provided, use default!',.TRUE.)
          RETURN
       ENDIF
       READ(WORDSOUT(NCURRENT),*,ERR=30) F
       RETURN
30     STOP 'input> Error while trying to read float'
    END SUBROUTINE READF

    SUBROUTINE READI(I)
       IMPLICIT NONE
       INTEGER :: I

       NCURRENT = NCURRENT + 1
       IF (NCURRENT.GT.NITEMS) THEN
          CALL REPORT('input> Too few arguments provided, use default!',.TRUE.)
          RETURN
       ENDIF
       READ(WORDSOUT(NCURRENT),*,ERR=40) I
       RETURN
40     STOP 'input> Error while trying to read integer'
    END SUBROUTINE READI

    SUBROUTINE READU(A)
       IMPLICIT NONE
       CHARACTER(*), INTENT(OUT) :: A
       
       NCURRENT = NCURRENT + 1
       IF (NCURRENT.GT.NITEMS) THEN
          CALL REPORT('input> Too few arguments provided, use default!',.TRUE.)
          RETURN
       ENDIF
       A=WORDSOUT(NCURRENT)
       CALL L2U(A)
       RETURN
    END SUBROUTINE READU

!Turns string into all upper case
    SUBROUTINE L2U(KEYWORD)
       IMPLICIT NONE
       CHARACTER(*),INTENT(INOUT) :: KEYWORD
       INTEGER :: LENGTH, I, IDX
       CHARACTER(LEN=26), PARAMETER :: U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       CHARACTER(LEN=26), PARAMETER :: L = 'abcdefghijklmnopqrstuvwxyz'     

       LENGTH=LEN(KEYWORD)

       DO I=1,LENGTH
          IDX = INDEX(L,KEYWORD(I:I))
          IF (IDX.GT.0) KEYWORD(I:I) = U(IDX:IDX)
       ENDDO
    END SUBROUTINE L2U  

!Turns string into all lower case
    SUBROUTINE U2L(KEYWORD)
       IMPLICIT NONE
       CHARACTER(*),INTENT(INOUT) :: KEYWORD
       INTEGER :: LENGTH, I, IDX
       CHARACTER(LEN=26), PARAMETER :: U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
       CHARACTER(LEN=26), PARAMETER :: L = 'abcdefghijklmnopqrstuvwxyz'     

       LENGTH=LEN(KEYWORD)

       DO I=1,LENGTH
          IDX = INDEX(U,KEYWORD(I:I))
          IF (IDX.GT.0) KEYWORD(I:I) = L(IDX:IDX)
       ENDDO
    END SUBROUTINE U2L 

    SUBROUTINE REPORT(C,REFLCT)
       IMPLICIT NONE
       CHARACTER(*) :: C
       LOGICAL      :: REFLCT
      
       WRITE(*,'(A)') C
       IF (REFLCT) THEN
          WRITE(*,'(A)') 'input> Current input line:'
          WRITE(*,'(A)') LINE         
       ENDIF
    END SUBROUTINE REPORT
 
    SUBROUTINE READ_LINE(NWORDS)
       IMPLICIT NONE
       INTEGER, INTENT(OUT) :: NWORDS

       CHARACTER(25), ALLOCATABLE :: TEMPOUT(:) 
       INTEGER:: J1,START_IND,END_IND,J2
       CHARACTER(25) :: WORD

       !remove content of previously read line
       IF (ALLOCATED(WORDSOUT)) DEALLOCATE(WORDSOUT)
       !reset all counter variables
       START_IND=0
       END_IND=0
       J1=1
       J2=0
       !Go through the line character by character
       DO WHILE(J1.LE.LEN(LINE))
          !while we don't know the start of the word, just read until we find a non-blank character
          IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
             START_IND=J1
          ENDIF
          !Found some charcters now - find the end of the string
          IF (START_IND.GT.0) THEN
             IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
             IF (J1.EQ.LEN(LINE)) END_IND=J1
             IF (END_IND.GT.0) THEN
                J2=J2+1
                !if it is the first word we have to allocate WORDSOUT
                IF (.NOT.ALLOCATED(WORDSOUT)) THEN
                   ALLOCATE(WORDSOUT(J2))
                !if it is not the first word we have to reallocate
                ELSE
                   ALLOCATE(TEMPOUT(J2-1))
                   TEMPOUT(:) = WORDSOUT(:)
                   DEALLOCATE(WORDSOUT)
                   ALLOCATE(WORDSOUT(J2))
                   WORDSOUT(1:J2-1) = TEMPOUT(1:J2-1)
                   DEALLOCATE(TEMPOUT)
                ENDIF
                WORD=LINE(START_IND:END_IND)
                WORDSOUT(J2)=TRIM(WORD)
                START_IND=0
                END_IND=0
             ENDIF
          ENDIF
          J1=J1+1
       ENDDO
       NWORDS=J2
    END SUBROUTINE READ_LINE

END MODULE INPUTmod
