MODULE MD_CALCS
   USE NUMKIND
   CONTAINS
      ! function to get kinetic energy
      REAL(KIND=REAL64) FUNCTION E_KINETIC(VEL,TEMP) RESULT(EKIN)
         USE MD_COMMONS, ONLY: NOPT, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: TEMP
         INTEGER :: I
         EKIN = 0.0D0
         DO I=1,NOPT
            EKIN = EKIN + VEL(I)*VEL(I)/MASSES(I)
         END DO
         EKIN = 0.5*EKIN
      END FUNCTION E_KINETIC

      ! subroutine to get com and linear momentum
      SUBROUTINE DETERMINE_LINMOM(X,VEL,COM,PCOM)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         REAL(KIND=REAL64), INTENT(OUT) :: COM(3)
         REAL(KIND=REAL64), INTENT(OUT) :: PCOM(3)
         INTEGER :: I, J, IDX
         REAL(KIND=REAL64) :: TOTALMASS

         COM(1:3) = 0.0D0
         PCOM(1:3) = 0.0D0
         DO I=1,NATOM
            DO J=1,3
               IDX = 3*(I-1) + J
               COM(J) = COM(J) + COORDS(IDX)*MASSES(J)
               P_COM(J) = P_COM(J) + VEL(IDX)*MASSES(J)              
            END DO
         END DO
         TOTALMASS = SUM(MASSES)/3
         COM = COM/TOTALMASS
      END SUBROUTINE DETERMINE_LINMOM

      SUBROUTINE REMOVE_LINMOM(X,VEL,CENTRET)
         USE MD_COMMONS, ONLY: NOPT, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         LOGICAL, INTENT(IN) :: CENTRET                 
         REAL(KIND=REAL64)

      END SUBROUTINE REMOVE_LINMOM
   ! 2. remove net translation and get rotation
   DO I=1,NATOMS
      R(1:3) = 0.0D0
      DO J=1,3 
         IDX = 3*(I-1) + J 
         ! change velocity to account for linear momentum
         VEL(IDX) = VEL(IDX) - P_COM(J)/TOTALMASS
         ! get position vector
         R(J) = COORDS(IDX) - COM(J)
         X(IDX) = R(J)
         ! get momentum vector
         P(J) = VEL(IDX) * MASSES(I)
      END DO
      R2 = DOT_PRODUCT(R,R)
      ANG_MOM = ANG_MOM + CROSSP(R, P)
      ! update moment of inertia
      DO J=1,3
         MOI(J,J) = MOI(J,J) + (R2 - R(J)*R(J))*MASSES(I)
         DO K=J,3
            MOI(J,K) = MOI(J,K) - R(J)*R(K)*MASSES(I)
            MOI(K,J) = MOI(J,K)
         END DO
      END DO
   END DO
END MODULE MD_CALCS