SUBROUTINE DUMPSTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)

USE COMMONS
USE PREC
USE PORFUNCS
USE OUTPUT, ONLY : WRITE_MARKOV_COORDS
IMPLICIT NONE
INTEGER NDONE, JBEST(NPAR), JP, J1, MYUNIT2, GETUNIT, LUNIT
REAL(KIND = REAL64) EBEST(NPAR), BESTCOORDS(3*NATOMSALLOC,NPAR)
CHARACTER(LEN=20) :: ISTR
CHARACTER(LEN=20) :: ISTR1
LOGICAL EXISTS

  IF (NPAR.GT.1) THEN
     WRITE (ISTR, '(i10)') JP
     MYUNIT2=GETUNIT()
     CALL FLUSH(MYUNIT)
     OPEN(MYUNIT2,FILE="GMIN.dump."//trim(adjustl(istr)),STATUS='UNKNOWN')
  ELSE
     MYUNIT2=GETUNIT()
     WRITE (ISTR, '(A)') ''
     INQUIRE (FILE='GMIN.dump', EXIST=EXISTS)
     IF (EXISTS) CALL SYSTEM('cp GMIN.dump GMIN.dump.save')
     OPEN(MYUNIT2,FILE='GMIN.dump',STATUS='UNKNOWN')
  ENDIF

  WRITE(MYUNIT2, '(A)') 'steps completed NQ(JP) in mc'
  WRITE(MYUNIT2, '(I8)') NDONE
  WRITE(MYUNIT2, '(I8)') NSAVE
! Number of potential energy calls done so far
  WRITE(MYUNIT2, '(I20)') NPCALL 
  WRITE(MYUNIT2, '(A)') 'COORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  WRITE(MYUNIT2, '(I8)') NATOMS
  CALL WRITE_MARKOV_COORDS(MYUNIT2, '(3F25.15)', JP)
!  WRITE(MYUNIT2, '(3F25.15)') COORDSO(1:3*NATOMS,JP)
  WRITE(MYUNIT2, '(A)') 'STEP, ASTEP, TEMP:'
  WRITE(MYUNIT2, '(3F25.15)') STEP(JP),ASTEP(JP),TEMP(JP)
  WRITE(MYUNIT2, '(A)') 'QMIN and QMINP'
! Dump info for each saved minimum
  DO J1=1,NSAVE
     WRITE(MYUNIT2, '(A,I8)') 'saved minimum ',J1
! Energy of saved minimum J1
     WRITE(MYUNIT2, '(G25.15)') QMIN(J1)
! Number of potential energy calls taken when this minimum was first encountered
     WRITE(MYUNIT2, '(I20)') NPCALL_QMIN(J1)
! Coordinates of saved minimum J1
     WRITE(MYUNIT2, '(I20)') QMINNATOMS(J1)
     WRITE(MYUNIT2, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1))
     IF (QMIN(J1).LT.1.0D10) THEN   
        WRITE (ISTR1, '(i10)') J1 
        LUNIT=GETUNIT()
        IF ((NPAR.GT.1).AND.(DUMPMINT)) THEN
           OPEN(LUNIT,FILE="dumpmin."//TRIM(ADJUSTL(ISTR))//"."//TRIM(ADJUSTL(ISTR1)),STATUS='UNKNOWN')
           WRITE(LUNIT, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1)) 
           CLOSE(LUNIT)
        ELSEIF (DUMPSTRUCTURES.AND.DUMPMINT) THEN 
           OPEN(LUNIT,FILE="dumpmin."//TRIM(ADJUSTL(ISTR1)),STATUS='UNKNOWN')
           WRITE(LUNIT, '(3F25.15)') QMINP(J1,1:3*QMINNATOMS(J1)) 
           CLOSE(LUNIT) 
        ENDIF
     ENDIF  
  ENDDO
  WRITE(MYUNIT2,'(A)') 'new restart procedure - JBEST, EBEST, BESTCOORDS'
  WRITE(MYUNIT2, '(A,I8)') 'run number ',JP
  WRITE(MYUNIT2, '(I8,F25.15)') JBEST(JP), EBEST(JP)
  WRITE(MYUNIT2, '(I8)') QMINNATOMS(1)
  WRITE(MYUNIT2, '(3F25.15)') BESTCOORDS(1:3*QMINNATOMS(1),JP)
  CLOSE(MYUNIT2)

RETURN
END SUBROUTINE DUMPSTATE

SUBROUTINE RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
USE COMMONS
USE PREC
IMPLICIT NONE
INTEGER NDONE, JBEST(NPAR), JP, J1, J2, OLDSAVE, MYUNIT2, GETUNIT, NDUMMY
REAL(KIND = REAL64) EBEST(NPAR), BESTCOORDS(3*NATOMS,NPAR)
CHARACTER(LEN=20) :: ISTR

IF (NPAR.GT.1) THEN
   WRITE (ISTR, '(i10)') JP
   MYUNIT2=GETUNIT()
   OPEN(MYUNIT2,FILE=TRIM(ADJUSTL(DUMPFILE))//trim(adjustl(istr)),STATUS='OLD')
ELSE
   MYUNIT2=1
   OPEN(UNIT=MYUNIT2,FILE=TRIM(ADJUSTL(DUMPFILE)),STATUS='OLD')
ENDIF

   READ(MYUNIT2,*) 
   READ(MYUNIT2,'(I8)') NDONE
   READ(MYUNIT2,'(I8)') OLDSAVE
   READ(MYUNIT2,'(I20)') NPCALL 
   READ(MYUNIT2,*) 
! Read in last minimum, which is also the last minimum in the Markov chain.
   READ(MYUNIT2,*) 
   READ(MYUNIT2,*) NATOMS
   READ(MYUNIT2,*) COORDS(1:3*NATOMS,JP)
   DO J1=1,NATOMS*(NPAR-1)
      READ(MYUNIT2,*)
   ENDDO

   WRITE(MYUNIT,'(A,I4)') 'Initial coordinates: process',JP
   WRITE(MYUNIT,'(3F20.10)') (COORDS(J1,JP),J1=1,3*NATOMS)

! step and temperature information
   READ(MYUNIT2,*) 
   READ(MYUNIT2,*) STEP(JP),ASTEP(JP),TEMP(JP)
! best OLDSAVE minima
   READ(MYUNIT2,*) 
!       this is the loop that has been edited to allow SAVE to vary between
!       runs. There are three cases which must be considered. 1) NSAVE=OLDSAVE.
!       In this case, nothing new needs to be done. 2) NSAVE<OLDSAVE. Now, the
!       extra minima that are saved in the dump file need to be cycled over
!       (ignored). 3) NSAVE>OLDSAVE. The array must be allocated with
!       dimensions suitable for the new number of minima, and the dump file read
!       in as far as possible. The rest of the array must then be padded with
!       zeros. The QMIN and QMINP arrays are already allocated in main.F with
!       the correct dimension so that need not worry us here :)
   IF (NSAVE.EQ.OLDSAVE) THEN
      DO J1=1,NSAVE
         READ(MYUNIT2,*) 
         READ(MYUNIT2,*) QMIN(J1)
         READ(MYUNIT2,*) NPCALL_QMIN(J1)
         READ(MYUNIT2,*) QMINNATOMS(J1)
         READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
      ENDDO
   ELSEIF (NSAVE.LT.OLDSAVE) THEN
      PRINT *,'NSAVE<OLDSAVE - truncating read in'
      DO J1=1,OLDSAVE
         IF (J1.LE.NSAVE) THEN
            READ(MYUNIT2,*) 
            READ(MYUNIT2,*) QMIN(J1)
            READ(MYUNIT2,*) NPCALL_QMIN(J1)
            READ(MYUNIT2,*) QMINNATOMS(J1)
            READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
         ELSE 
            
            READ(MYUNIT2,*) 
            READ(MYUNIT2,*) 
            READ(MYUNIT2,*) 
            READ(MYUNIT2,*) NDUMMY
            DO J2=1,NDUMMY
               READ(MYUNIT2,*)
            ENDDO
         ENDIF
      ENDDO
   ELSEIF (NSAVE.GT.OLDSAVE) THEN
      PRINT *,'NSAVE>OLDSAVE - padding QMIN and QMINP with zeros'
      DO J1=1,NSAVE
         IF (J1.LE.OLDSAVE) THEN
            READ(MYUNIT2,*) 
            READ(MYUNIT2,*) QMIN(J1)
            READ(MYUNIT2,*) NPCALL_QMIN(J1)
            READ(MYUNIT2,*) QMINNATOMS(J1)
            READ(MYUNIT2,*) QMINP(J1,1:3*QMINNATOMS(J1))
         ELSE 
            QMIN(J1)=0.0D0
            NPCALL_QMIN(J1)=0
            QMINNATOMS(J1)=NATOMS
            QMINP(J1,1:3*QMINNATOMS(J1))=0.0D0
         ENDIF
      ENDDO
   ENDIF
       
   ! read in EBEST and BESTCOORDS
   READ(MYUNIT2,*) 
   READ(MYUNIT2,*) 
   READ(MYUNIT2,*) JBEST(JP), EBEST(JP)
   READ(MYUNIT2,*) NDUMMY
   READ(MYUNIT2,*) BESTCOORDS(1:3*NDUMMY,JP)
   CLOSE(MYUNIT2)

RETURN
END SUBROUTINE RESTORESTATE
