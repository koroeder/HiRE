MODULE MODHESS
      USE PREC
      IMPLICIT NONE
      SAVE

      REAL(KIND = REAL64), ALLOCATABLE :: HESS(:, :)  !  3*MXATMS,3*MXATMS

      ! hk286 - toggle between moving and stationary frame
      LOGICAL :: RBAANORMALMODET
      LOGICAL :: MASS_WEIGHTED = .FALSE.

CONTAINS

SUBROUTINE MASSWT()
      USE COMMONS, ONLY: NATOMS, ATMASS

      IMPLICIT NONE

      INTEGER :: J1, J2
      INTEGER :: X1, Y1, Z1, X2, Y2, Z2
      REAL(KIND = REAL64) :: AMASS, BMASS, FMASS

! Only apply mass weighting if it hasn't already been applied.
! Currently disabled.
!      IF (.NOT. MASS_WEIGHTED) THEN
          DO J1 = 1, NATOMS
             X1 = 3 * J1 - 2
             Y1 = 3 * J1 - 1
             Z1 = 3 * J1
             AMASS = 1.D0/SQRT(ATMASS(J1))
             DO J2 = J1, NATOMS
                X2 = 3 * J2 - 2
                Y2 = 3 * J2 - 1
                Z2 = 3 * J2
                BMASS = 1.D0/SQRT(ATMASS(J2))
                FMASS = AMASS*BMASS
                HESS(X1,X2) = FMASS*HESS(X1,X2)
                HESS(X1,Y2) = FMASS*HESS(X1,Y2)
                HESS(X1,Z2) = FMASS*HESS(X1,Z2)
                HESS(Y1,X2) = FMASS*HESS(Y1,X2)
                HESS(Y1,Y2) = FMASS*HESS(Y1,Y2)
                HESS(Y1,Z2) = FMASS*HESS(Y1,Z2)
                HESS(Z1,X2) = FMASS*HESS(Z1,X2)
                HESS(Z1,Y2) = FMASS*HESS(Z1,Y2)
                HESS(Z1,Z2) = FMASS*HESS(Z1,Z2)
                IF (J1 .NE. J2) THEN
                   HESS(X2,X1) = HESS(X1,X2)
                   HESS(Y2,X1) = HESS(X1,Y2)
                   HESS(Z2,X1) = HESS(X1,Z2)
                   HESS(X2,Y1) = HESS(Y1,X2)
                   HESS(Y2,Y1) = HESS(Y1,Y2)
                   HESS(Z2,Y1) = HESS(Y1,Z2)
                   HESS(X2,Z1) = HESS(Z1,X2)
                   HESS(Y2,Z1) = HESS(Z1,Y2)
                   HESS(Z2,Z1) = HESS(Z1,Z2)
                ENDIF
             ENDDO
          ENDDO
!          MASS_WEIGHTED = .TRUE.
!      END IF    
END SUBROUTINE MASSWT

END MODULE MODHESS
