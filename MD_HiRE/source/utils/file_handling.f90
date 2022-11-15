!> @file
!> Contains FILE_UTILS module with utilities to handle files

!> Module to handle basic file operations
!> @brief
!> 
!> This module contains utility functions used across PiGS.
!> The functionality in the module includes opening files, getting empty units, 
!> obtaining the length of files, checking whether a file or directory exists, and 
!> printing the current directory.

MODULE FILE_UTILS
   ! File_open, Getunit and File_length are taken from GMIN 
   IMPLICIT NONE
   CONTAINS

      !> Subroutine to open file
      !> @brief
      !>
      !> The subroutine uses GETUNIT to find an empty unit and then opens the file on this unit.
      !> The file can be opened to append or not. 
      !> The subroutine does not check whether the file in question exists.
      !>
      !> @param[in] FILE_NAME - name of file to be opened (string)     
      !> @param[out] FILE_UNIT - unit used to open file (integer)
      !> @param[in] APPEND - Open file to be appended (logical)
      !>
      !> @return Unit associated with opended file
      !>
      !> @see GETUNIT  
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

      !> Function to find empty unit to open file
      !>
      !> @return Free file unit
      !>
      !> @warning The unit returned is free at the time of the call, but is not reserved.
      !> A GETUNIT call should be immediately followed by an OPEN call. Alternatively, use FILE_OPEN.
      !> 
      !> @see FILE_OPEN
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

      !> Function to measure length of file
      !> @brief
      !>
      !> The file provided is opened with FILE_OPEN.
      !> The file is then read line by line, and a counter increased for every line. 
      !> IOSTAT is used to validate whether the file is read completely.
      !>
      !> @param[in] FILE_NAME - file for which we want to know the file length
      !>
      !> @return Number of lines in file
      !> 
      !> @see FILE_OPEN
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
      
      !> Function to check whether a directory exists
      !> @brief
      !>
      !> The directory name is checked using a system call.
      !> We use "test -d" to check the existence.
      !>
      !> @param[in] DIRNAME - directory to be checked
      !>
      !> @return True/False for existence of directory
      LOGICAL FUNCTION DIR_EXIST(DIRNAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: DIRNAME
         INTEGER :: STAT
         CALL EXECUTE_COMMAND_LINE("test -d "//DIRNAME,EXITSTAT=STAT)
         IF (STAT.EQ.0) THEN
            DIR_EXIST = .TRUE.
         ELSE
            DIR_EXIST = .FALSE.
         ENDIF
      END FUNCTION DIR_EXIST
      
      !> Function to check whether a file exists
      !> @brief
      !>
      !> The file name is checked using INQUIRE.
      !>
      !> @param[in] FILENAME - directory to be checked
      !>
      !> @return True/False for existence of file
      LOGICAL FUNCTION FILE_EXIST(FILENAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: FILENAME 
         INQUIRE(FILE=FILENAME, EXIST=FILE_EXIST)
      END FUNCTION FILE_EXIST 
      
!      !> Prints current directory for Debugging
!      !>
!      !> @param[in] OUTUNIT - Unit used for output to be written
!      !>
!      !> @return Prints current directory in file attached to unit provided1
!      SUBROUTINE PRINT_CURRDIR(OUTUNIT)
!         IMPLICIT NONE
!         INTEGER, INTENT(IN) :: OUTUNIT
!         CHARACTER(LEN=255) :: CWD
!         
!         CALL GETCWD(CWD)
!         WRITE(OUTUNIT,'(2A)') " current directory: ", TRIM(ADJUSTL(CWD))   
!      END SUBROUTINE PRINT_CURRDIR 
      
END MODULE FILE_UTILS
