!   Writes the output file at the end of a GMIN run
!
SUBROUTINE FINALIO
  USE COMMONS
  USE PREC
  USE GENRIGID, ONLY : RIGIDINIT
  USE HIRE_INTERFACE, ONLY: TERMINATE_HIRE, DUMP_PDB, PASS_PARTICLE_NAMES

  IMPLICIT NONE

  INTEGER             :: J1, J2, GETUNIT, MYUNIT2, MYUNIT3
  REAL(KIND = REAL64) :: TEND
  CHARACTER(LEN=20)   :: ISTR, MYFILENAME2
  CHARACTER(LEN=6)    :: J1_STRING, MYNODE_STRING
  CHARACTER(LEN=4)    :: ATOM_NAMES(NATOMS)

  CALL PASS_PARTICLE_NAMES(NATOMS,ATOM_NAMES)

  IF (MPIT) THEN
     WRITE (ISTR, '(I10)') MYNODE+1
     MYUNIT2=GETUNIT()
     MYFILENAME2="lowest."//TRIM(ADJUSTL(ISTR))
     OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(MYFILENAME2)), STATUS="UNKNOWN", FORM="FORMATTED")
  ELSE
     MYUNIT2=GETUNIT()
     IF (RIGIDINIT) THEN
        OPEN(MYUNIT2,FILE='GRlowest',STATUS='UNKNOWN')
     ELSE
        OPEN(MYUNIT2,FILE='lowest',STATUS='UNKNOWN')
     ENDIF
  ENDIF

  savemin: DO J1=1,NSAVE
     WRITE(MYUNIT2,'(I8)') NATOMS
     WRITE(MYUNIT2,10) J1, QMIN(J1), FF(J1), NPCALL_QMIN(J1)
10   FORMAT('Energy of minimum ',I6,'=',G20.10,' first found at step ',I8,' after ',I20,' function calls')
     DO J2=1,NATOMS
        WRITE(MYUNIT2,'(A4,3F16.7)') ATOM_NAMES(J2),QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
     ENDDO
     !we only have pdb files to write
     WRITE(J1_STRING,'(I6)') J1
     IF (MPIT) THEN
        WRITE(MYNODE_STRING,'(I6)') MYNODE + 1
        CALL DUMP_PDB(3*NATOMS, QMINP(J1,:),'lowest.'//TRIM(ADJUSTL(J1_STRING))//&
                   &'.'//TRIM(ADJUSTL(MYNODE_STRING))//'.pdb',.TRUE.)
        !also write start file to be used for OPTIM runs
        MYUNIT3=GETUNIT()
        MYFILENAME2='start.'//TRIM(ADJUSTL(J1_STRING))//'.'//trim(adjustl(MYNODE_STRING))
        OPEN(MYUNIT3,FILE=trim(adjustl(MYFILENAME2)), STATUS="unknown",form="formatted")
        DO J2=1,NATOMS
           WRITE(MYUNIT3,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
        ENDDO
        CLOSE(MYUNIT3)
     ELSE
        CALL DUMP_PDB(3*NATOMS, QMINP(J1,:),'lowest.'//TRIM(ADJUSTL(J1_STRING))//'.pdb',.TRUE.)
        !also write start file to be used for OPTIM runs
        MYUNIT3=GETUNIT()
        MYFILENAME2='start.'//TRIM(ADJUSTL(J1_STRING))
        OPEN(MYUNIT3,FILE=trim(adjustl(MYFILENAME2)),STATUS="unknown",form="formatted")
        DO J2=1,NATOMS
           WRITE(MYUNIT3,'(3F28.20)') QMINP(J1,3*(J2-1)+1),QMINP(J1,3*(J2-1)+2),QMINP(J1,3*(J2-1)+3)
        ENDDO
        CLOSE(MYUNIT3)
     ENDIF
     
  ENDDO savemin

  CLOSE(MYUNIT2)

  CALL TERMINATE_HIRE()

  CALL CPU_TIME(TEND)
  WRITE(MYUNIT,"(A,F18.1,A)") "time elapsed ", TEND - TSTART, " seconds"
  WRITE(MYUNIT,"(A,I18)") "Number of potential calls ", NPCALL
  RETURN

END SUBROUTINE FINALIO


