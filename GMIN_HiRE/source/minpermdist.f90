!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  COORDSA becomes the optimal alignment of the optimal permutation(-inversion)
!  isomer, but without the permutations. DISTANCE is the residual square distance
!  for the best alignment with respect to permutation(-inversion)s as well as
!  orientation and centre of mass.


SUBROUTINE MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,DISTANCE,RMATBEST)
USE PREC
USE GENRIGID
USE PORFUNCS
IMPLICIT NONE

!Arguments
INTEGER             :: NATOMS
REAL(KIND = REAL64) :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, RMATBEST(3,3) 
LOGICAL             :: DEBUG

!Variables
REAL(KIND = REAL64) :: RMAT(3,3), TEMPCOORDSA(DEGFREEDOMS), TEMPCOORDSB(DEGFREEDOMS)

DISTANCE = 0.0D0
RMAT(:,:) = 0.0D0

IF (RIGIDINIT) THEN
    IF(DEBUG) THEN
        IF(.NOT.(ANY(ABS(COORDSA(DEGFREEDOMS+1:3*NATOMS)) .GT. 1.0E-10))) THEN
            WRITE(*,*) "minpermdist> Warning: COORDSA seems to be in AA coords. Last block (should all be 0):"
            WRITE(*,*) COORDSA(DEGFREEDOMS+1:3*NATOMS)
            WRITE(*,*) "Transforming to Cartesians."
            TEMPCOORDSA = COORDSA(:DEGFREEDOMS)
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, TEMPCOORDSA, COORDSA)
            TEMPCOORDSA(:) = 0
        ENDIF
        IF(.NOT.(ANY(ABS(COORDSB(DEGFREEDOMS+1:3*NATOMS)) .GT. 1.0E-10))) THEN
            WRITE(*,*) "minpermdist> Warning: COORDSB seems to be in AA coords. Last block (should all be 0):"
            WRITE(*,*) COORDSB(DEGFREEDOMS+1:3*NATOMS)
            WRITE(*,*) "Transforming to Cartesians."
            TEMPCOORDSB = COORDSB(:DEGFREEDOMS)
            CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, TEMPCOORDSB, COORDSB)
            TEMPCOORDSB(:) = 0
        ENDIF
    ENDIF
ENDIF

CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,DEBUG,RMAT)
RMATBEST = RMAT
RETURN
END SUBROUTINE MINPERMDIST
