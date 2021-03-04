      SUBROUTINE MYINERTIA(COORDS,NATOMS,ITDET,I1,I2,I3)
      USE COMMONS, ONLY: ATMASS
      USE PREC
      IMPLICIT NONE
      INTEGER J1, J2, J3, NATOMS
      REAL(KIND = REAL64) ITDET,COORDS(3*NATOMS),DUMQ(3*NATOMS),I1,I2,I3
      REAL(KIND = REAL64) IT(3,3), CMX, CMY, CMZ, MASST, IV(3,3)

      DUMQ(1:3*NATOMS)=COORDS(1:3*NATOMS)
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      MASST=0.0D0
      DO J1=1,NATOMS
            CMX=CMX+DUMQ(3*(J1-1)+1)*ATMASS(J1)
            CMY=CMY+DUMQ(3*(J1-1)+2)*ATMASS(J1)
            CMZ=CMZ+DUMQ(3*(J1-1)+3)*ATMASS(J1)
            MASST=MASST+ATMASS(J1)
      ENDDO
      CMX=CMX/MASST
      CMY=CMY/MASST
      CMZ=CMZ/MASST
      DO J1=1,NATOMS
         DUMQ(3*(J1-1)+1)=DUMQ(3*(J1-1)+1)-CMX
         DUMQ(3*(J1-1)+2)=DUMQ(3*(J1-1)+2)-CMY
         DUMQ(3*(J1-1)+3)=DUMQ(3*(J1-1)+3)-CMZ
      ENDDO

      DO J1=1,3
         DO J2=1,3
            IT(J1,J2)=0.0D0
            j3loop1: DO J3=1,NATOMS
               IT(J1,J2)=IT(J1,J2)-DUMQ(3*(J3-1)+J1)*DUMQ(3*(J3-1)+J2)*ATMASS(J3)
            ENDDO j3loop1
            IF (J1.EQ.J2) THEN
               j3loop2: DO J3=1,NATOMS
                  IT(J1,J2)=IT(J1,J2)+(DUMQ(3*(J3-1)+1)**2+DUMQ(3*(J3-1)+2)**2+DUMQ(3*(J3-1)+3)**2)*ATMASS(J3)
               ENDDO j3loop2
            ENDIF
         ENDDO
      ENDDO

! Diagonalize inertia tensor. The 0 flag reorders the e/values and
! e/vectors from smallest to largest.
      CALL EIG(IT,IV,3,3,0)
      I1=IT(1,1)
      I2=IT(2,2)
      I3=IT(3,3)

      ITDET=IT(1,1)*IT(2,2)*IT(3,3)
      IF (NATOMS.EQ.2) ITDET=IT(2,2)*IT(3,3)

      RETURN
      END

