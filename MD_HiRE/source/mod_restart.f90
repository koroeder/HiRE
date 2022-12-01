MODULE MOD_RESTART
   USE NUMKIND
   IMPLICIT NONE
   CHARACTER(LEN=25) :: RSTNAME = "md_restart.dat"
   CHARACTER(LEN=25) :: RSTNAMEOLD = "md_restart_old.dat"   
   CONTAINS
      SUBROUTINE WRITE_RST_FILE(CURRSTEP, X, VEL)
         USE MD_COMMONS, ONLY: NATOMS, NOPT, TEMP
         USE FILE_UTILS, ONLY: FILE_EXIST, FILE_OPEN
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND = REAL64), INTENT(IN) :: VEL(NOPT)
         INTEGER :: RSTUNIT

         IF (FILE_EXIST(RSTNAME)) THEN
            CALL EXECUTE_COMMAND_LINE("mv "//TRIM(ADJUSTL(RSTNAME))//" "//TRIM(ADJUSTL(RSTNAMEOLD)))
         END IF
         CALL FILE_OPEN(RSTNAME, RSTUNIT, .FALSE.)
         WRITE(RSTUNIT,'(A,I8)') "NATOMS ", NATOMS
         WRITE(RSTUNIT,'(A,I10)') "STEP ", CURRSTEP
         WRITE(RSTUNIT,'(A,F9.4)') "TEMP ", TEMP         
         WRITE(RSTUNIT,'(3(F15.7))') X
         WRITE(RSTUNIT,'(3(F15.7))') VEL
         CLOSE(RSTUNIT)               
      END SUBROUTINE WRITE_RST_FILE

      SUBROUTINE READ_RST_FILE(INPF, CURRSTEP, TEMP, X, VEL)
         USE MD_COMMONS, ONLY: MYUNIT, NATOMS, NOPT
         USE FILE_UTILS, ONLY: FILE_EXIST, FILE_OPEN
         USE INPUTMOD, ONLY: INPUTKW, READI, READF
         USE MD_UTILS, ONLY: TERMINATE_ERR
         IMPLICIT NONE 
         CHARACTER(LEN=25), INTENT(IN) :: INPF                      
         INTEGER, INTENT(OUT) :: CURRSTEP
         REAL(KIND = REAL64), INTENT(OUT) :: X(NOPT)
         REAL(KIND = REAL64), INTENT(OUT) :: VEL(NOPT)         
         REAL(KIND = REAL64), INTENT(OUT) :: TEMP
         INTEGER :: RSTUNIT, NATS, J
         CHARACTER(LEN=20) :: KEYWORD
         LOGICAL :: ENDT

         IF (.NOT.FILE_EXIST(INPF)) THEN
            WRITE(MYUNIT,*) " read_rst> File ", INPF, " does not exist - STOP"
            CALL TERMINATE_ERR(.TRUE.,.TRUE.)
         END IF

         CALL FILE_OPEN(INPF,RSTUNIT,.FALSE.)
         ! first line is NATOMS
         CALL INPUTKW(ENDT,KEYWORD,RSTUNIT,.FALSE.)
         IF (ENDT) THEN
            CLOSE(RSTUNIT)
            WRITE(MYUNIT,'(A)') " readrst> Restart file incomplete - STOP"
            CALL TERMINATE_ERR(.TRUE.,.TRUE.)
         ELSE
            CALL READI(NATS)
            IF (NATS.NE.NATOMS) THEN
                CLOSE(RSTUNIT)
                WRITE(MYUNIT,'(A)') " readrst> Number of atoms in restart file is different from NATOMs from potential: ", &
                                    NATS, NATOMS,  "- STOP"
                CALL TERMINATE_ERR(.TRUE.,.TRUE.)
            END IF
         END IF
         ! second line is number of steps
         CALL INPUTKW(ENDT,KEYWORD,RSTUNIT,.FALSE.)
         IF (ENDT) THEN
            CLOSE(RSTUNIT)
            WRITE(MYUNIT,'(A)') " readrst> Restart file incomplete - STOP"
            CALL TERMINATE_ERR(.TRUE.,.TRUE.)
         ELSE
            CALL READI(CURRSTEP)
         END IF
         ! third line is temperature
         CALL INPUTKW(ENDT,KEYWORD,RSTUNIT,.FALSE.)
         IF (ENDT) THEN
            CLOSE(RSTUNIT)
            WRITE(MYUNIT,'(A)') " readrst> Restart file incomplete - STOP"
            CALL TERMINATE_ERR(.TRUE.,.TRUE.)
         ELSE
            CALL READF(TEMP)
         END IF
         ! now the coordinates
         READ(RSTUNIT,'(3F15.7)') (X(J), J=1,NOPT)
         ! finally the velocity
         READ(RSTUNIT,'(3F15.7)') (VEL(J), J=1,NOPT)         
         CLOSE(RSTUNIT)
      END SUBROUTINE READ_RST_FILE   

END MODULE MOD_RESTART