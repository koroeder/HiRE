!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.

      SUBROUTINE GSORT2()
      USE COMMONS
      USE PREC
      IMPLICIT NONE
      INTEGER J1, J2, L, NTEMP
      REAL(KIND = REAL64) DUMMY, C

      DO 20 J1=1,NSAVE-1
         L=J1
         DO 10 J2=J1+1,NSAVE
            IF (QMIN(L).GT.QMIN(J2)) L=J2
10       CONTINUE
         DUMMY=QMIN(L)
         QMIN(L)=QMIN(J1)
         QMIN(J1)=DUMMY
         NTEMP=QMINNATOMS(L)
         QMINNATOMS(L)=QMINNATOMS(J1)
         QMINNATOMS(J1)=NTEMP
! Swap the first found (FF) and number of potential calls when first found (NPCALL_QMIN) array elements to match
         NTEMP=FF(L)
         FF(L)=FF(J1)
         FF(J1)=NTEMP
         NTEMP=NPCALL_QMIN(L)
         NPCALL_QMIN(L)=NPCALL_QMIN(J1)
         NPCALL_QMIN(J1)=NTEMP
! Swap the coordinates to match
         DO J2=1,3*NATOMS
            C=QMINP(L,J2)
            QMINP(L,J2)=QMINP(J1,J2)
            QMINP(J1,J2)=C
         ENDDO
20    CONTINUE

      IF (FEBHT) THEN
         DO J1 = NSAVE, 1, -1
            QMIN(J1) = QMIN(J1) - QMIN(1)
         ENDDO
      ENDIF

      RETURN
      END
