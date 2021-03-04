!   This subroutine will be used to do all initializations, e.g. prepare
!   variables/arrays used in potentials, take step routines, etc.  It is assumed
!   that at this point the keywords have been read and the coords have been
!   loaded.

SUBROUTINE INITIALIZATIONS()
USE COMMONS
IMPLICIT NONE
INTEGER J1, J2, J6

! set HARMONICR0 to be the initial coords
IF ( HARMONICF ) THEN
   HARMONICR0(1:3*NATOMS) = COORDS(1:3*NATOMS,1)
ENDIF

! create FROZENLIST from FROZEN
! FROZENLIST holds the sorted list of frozen atoms between 1 and NFREEZE.
! Between NFREEZE+1 and N it holds the sorted list of unfrozen particles.
IF ( FREEZE ) THEN
   ALLOCATE(FROZENLIST(NATOMS))
   J1=0
   J2=NFREEZE
   DO J6=1,NATOMS
      IF ( FROZEN(J6) ) THEN
         J1 = J1+1
         FROZENLIST(J1) = J6
      ELSE
         J2 = J2+1
         FROZENLIST(J2) = J6
      ENDIF
   END DO
ENDIF



RETURN
END SUBROUTINE INITIALIZATIONS
