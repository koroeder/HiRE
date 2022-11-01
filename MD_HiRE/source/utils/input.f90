!> @file
!> Contains INPUTMOD module with utilities to read lines and parse them from opened files

!> Module to handle input
!> @brief
!>
!> Originally implemented for GMIN by K. Roeder \n 
!> Dependency on NUMKIND \n 
!> Subroutines available: \n 
!> INPUT  \n 
!>  -> returns keyword and detects end of file \n 
!>  -> input required: file unit and whether capatilisation is required \n 
!>  -> sets the internal variables NITEMS and the list WORDSOUT \n 
!>  -> calls READ_LINE and L2U \n 
!> READA, READI, READF and READU \n 
!>  -> return item in WORDSOUT as string, int, float, capatilised string \n 
!>  -> no input required, returns value into correct format \n 
!>  -> if there is no more item in the list, the subroutines return a warning
!>     and the default value for the variable is retained \n 
!>  -> requires that INPUT was called before(!), calls to REPORT and L2U \n 
!> Internal subroutines: \n 
!>  -> L2U and U2L: lower to upper and upper to lower case transformation for strings \n 
!>  -> REPORT: report errors with additional detail \n 
!>  -> READ_LINE: read line and return list of entries, respects " " and ' ' 
MODULE INPUTmod

   USE NUMKIND

   PRIVATE
   PUBLIC :: INPUT, READA, READF, READI, READU, NITEMS, INPUTKW, REPORT, RETURNENTRY
   !> number of items in current input line
   INTEGER :: NITEMS   
   !> current item read
   INTEGER :: NCURRENT  
   !> current input line from file
   CHARACTER (LEN=200) :: LINE
   !>list of all inputs (all items are strings)
   CHARACTER(80), ALLOCATABLE :: WORDSOUT(:)  
   !>Set debug mode for input module, this should be used for testing only  
   LOGICAL :: DEBUGMOD = .FALSE.

   CONTAINS

      !> Subroutine to read input line that contains a keyword
      !> @brief
      !>
      !> The first step is that we read a new line from the file unit provided.
      !> We then call READ_LINE to parse the line into a list of entries.
      !> If there are no entries, we take the next line.
      !> Otherwise, we return the first item as keyword and set the module variables accordingly.
      !> If UPT is true, we call L2U on keyword.
      !>
      !> @param[in] IR - file unit to read from     
      !> @param[in] UPT - capitalise the keyword?
      !> @param[out] ENDT - end of file status
      !> @param[out] KEYWORD - keyword found at start of input line
      !>
      !> @return Keyword and end of file status
      !>
      !> @see READ_LINE
      !> @see L2U
      SUBROUTINE INPUTKW(ENDT,KEYWORD,IR,UPT)
         IMPLICIT NONE
         LOGICAL, INTENT(OUT)      :: ENDT     !end of file?
         CHARACTER(*), INTENT(OUT) :: KEYWORD  !first word in line
         INTEGER, INTENT(IN)       :: IR       !file unit
         LOGICAL, INTENT(IN)       :: UPT      !logical to capatilise KEYWORD 

         INTEGER :: NENTRIES, LENGTH, I
         
10       READ(IR,'(A)',END=20,ERR=25) LINE
         NITEMS = 0
         CALL READ_LINE(NENTRIES)
         IF (NENTRIES.EQ.0) THEN
            IF (DEBUGMOD) WRITE(*,'(A)') 'input> Blank line in data file, reading next line'
            GOTO 10
         ELSE
            NITEMS = NENTRIES
            KEYWORD = WORDSOUT(1)
            IF (DEBUGMOD) WRITE(*,'(3A,I8)') ' input> Number of options specified for keyword ' &
                                             , KEYWORD, ' is ', NITEMS
         ENDIF
         IF (UPT) CALL L2U(KEYWORD)
         ENDT=.FALSE.
         NCURRENT=1
         RETURN
20       ENDT=.TRUE.
         RETURN   
25       STOP 'input> Error while reading next line'
      END SUBROUTINE


      !> Subroutine to read input line that contains no keyword
      !> @brief
      !>
      !> This routine works like INPUT_KW, but the first entry is not returned as keyword.
      !> Input information is stored in the module variables.
      !>
      !> @param[in] IR - file unit to read from     
      !> @param[out] ENDT - end of file status
      !>
      !> @return End of file status
      !>
      !> @see READ_LINE     
      SUBROUTINE INPUT(ENDT,IR)
         IMPLICIT NONE
         LOGICAL, INTENT(OUT)      :: ENDT     !end of file?
         INTEGER, INTENT(IN)       :: IR       !file unit

         INTEGER :: NENTRIES, LENGTH, I
         
10       READ(IR,'(A)',END=20,ERR=25) LINE
         NITEMS = 0
         CALL READ_LINE(NENTRIES)
         IF (NENTRIES.EQ.0) THEN
            IF (DEBUGMOD) WRITE(*,'(A)') 'input> Blank line in file, reading next line'
            GOTO 10
         ELSE
            NITEMS = NENTRIES
            IF (DEBUGMOD) WRITE(*,'(A,I8)') ' input> Number of entries is ', NITEMS
         ENDIF
         ENDT=.FALSE.
         NCURRENT=0
         RETURN
20       ENDT=.TRUE.
         RETURN   
25       STOP 'input> Error while reading next line'
      END SUBROUTINE

      !> Subroutine to read next item as string
      !> @brief
      !>
      !> Calls the next item from WORDSOUT, which was set by INPUT or INPUT_KW.
      !> Uses REPORT if there not enough items, and the return is the default value of the variable.
      !>
      !> @param[out] A - Read item form list as string   
      !>
      !> @return Next item from input line
      !>
      !> @warning If there not enough items, this is reported, but the program does not abort.
      !>
      !> @see REPORT     
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


      !> Subroutine to read next item as float
      !> @brief
      !>
      !> Calls the next item from WORDSOUT, which was set by INPUT or INPUT_KW.
      !> Uses REPORT if there not enough items, and the return is the default value of the variable.
      !>
      !> @param[out] F - Read item form list as float   
      !>
      !> @return Next item from input line
      !>
      !> @warning If there not enough items, this is reported, but the program does not abort.
      !>
      !> @see REPORT    
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
30       STOP 'input> Error while trying to read float'
      END SUBROUTINE READF


      !> Subroutine to read next item as integer
      !> @brief
      !>
      !> Calls the next item from WORDSOUT, which was set by INPUT or INPUT_KW.
      !> Uses REPORT if there not enough items, and the return is the default value of the variable.
      !>
      !> @param[out] I - Read item form list as integer  
      !>
      !> @return Next item from input line
      !>
      !> @warning If there not enough items, this is reported, but the program does not abort.
      !>
      !> @see REPORT    
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
40       STOP 'input> Error while trying to read integer'
      END SUBROUTINE READI


      !> Subroutine to read next item as uppercase string
      !> @brief
      !>
      !> Calls the next item from WORDSOUT, which was set by INPUT or INPUT_KW.
      !> Uses REPORT if there not enough items, and the return is the default value of the variable.
      !>
      !> @param[out] U - Read item form list as string   
      !>
      !> @return Next item from input line
      !>
      !> @warning If there not enough items, this is reported, but the program does not abort.
      !>
      !> @see REPORT    
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

      !> Subroutine to turn lower case into upper case letters
      !> @brief
      !>
      !> Acts on string to turn all lower case letters into upper case.
      !> The input string is manipulated and only lower case letters are altered.
      !> INDEX is used to find all lower case letters and replaced.
      !>
      !> @param KEYWORD - String to be manipulated  
      !>
      !> @return Upper case string
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

      !> Subroutine to turn upper case into lower letters
      !> @brief
      !>
      !> Acts on string to turn all upper case letters into lower case.
      !> The input string is manipulated and only upper case letters are altered.
      !> INDEX is used to find all upper case letters and replaced.
      !>
      !> @param KEYWORD - String to be manipulated  
      !>
      !> @return Lower case string
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

      !> Reporter function
      !> @brief
      !>
      !> Subroutine to reflect on input reading issues
      !> The item taht caused the issue is printed, and if REFLCT is turned on the current input line is also provided.
      !>
      !> @param[out] C - Item that we report on
      !> @param[in] REFLCT - Pritn additional information?    
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
   
      !> Subroutine to return entry by index
      !> @brief
      !>
      !> Allows us access to an item in the WORDSOUT list by index.
      !> The READX routines all go in order, and this provide an option to access any item in the input line.
      !> 
      !> @param[in] IENTRY - Index of item we want to obtain
      !> @param[out] STRING - Item returned as string     
      SUBROUTINE RETURNENTRY(IENTRY,STRING)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IENTRY 
         CHARACTER(LEN=80) :: STRING
         
         STRING = WORDSOUT(IENTRY)

      END SUBROUTINE RETURNENTRY
   
      !> Read line from input file
      !> @brief
      !>
      !> Key function to handle inpuit in this module.
      !> The line saved in the module variables is taken and processed into a list of substrings.
      !> We use a list that we reallocate as we go along.
      !> Every time we have a new word to save we reallocate WORDSOUT.
      !> The number of words is coutned, and this is what is returned.
      !> The subroutine accepts single and double quotes.
      !> We separate words by spaces, but single and double quotes override this and act as highger priority delimiters.
      !> 
      !> @param[out] NWORDS - Number of words in the input line      
      SUBROUTINE READ_LINE(NWORDS)
         IMPLICIT NONE
         INTEGER, INTENT(OUT) :: NWORDS
      
         LOGICAL :: SQUOTET, DQUOTET
         CHARACTER(80), ALLOCATABLE :: TEMPOUT(:) 
         INTEGER:: J1,START_IND,END_IND,J2
         CHARACTER(80) :: WORD

         SQUOTET = .FALSE.
         DQUOTET = .FALSE.
   
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
               !if the non-blank character is a quote we want to read up to the next quote, otherwise we aim for the next blank
               !if we have quotes, the first character is the next one, not the current!
               IF (LINE(J1:J1).EQ.'"') THEN
                  DQUOTET=.TRUE.
                  START_IND=J1+1
               ELSE IF (LINE(J1:J1).EQ."'") THEN
                  SQUOTET=.TRUE.
                  START_IND=J1+1
               ELSE
                  START_IND=J1
               ENDIF
            !Found some characters now - find the end of the string
            ELSE IF (START_IND.GT.0) THEN
               IF (DQUOTET) THEN
                  IF (LINE(J1:J1).EQ.'"') END_IND=J1-1
               ELSE IF (SQUOTET) THEN
                  IF (LINE(J1:J1).EQ."'") END_IND=J1-1
               ELSE
                  IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
               ENDIF
               IF (J1.EQ.LEN(LINE)) END_IND=J1
               IF ((DQUOTET.AND.SQUOTET).AND.(START_IND.NE.0)) THEN
                  !quote followed by quote (empty) - clearing and move on
                  IF (START_IND.EQ.END_IND) THEN
                     DQUOTET=.FALSE.
                     SQUOTET=.FALSE.
                     START_IND=0
                     END_IND=0
                  ENDIF
               ENDIF
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
                  DQUOTET=.FALSE.
                  SQUOTET=.FALSE.
                  START_IND=0
                  END_IND=0
               ENDIF
            ENDIF
            J1=J1+1
         ENDDO
         NWORDS=J2
      END SUBROUTINE READ_LINE

END MODULE INPUTmod
