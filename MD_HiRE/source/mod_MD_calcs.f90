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
      SUBROUTINE DETERMINE_LINMOM(VEL)

      END SUBROUTINE DETERMINE_LINMOM

      SUBROUTINE REMOVE_LINMOM(X,VEL,CENTRET)
         USE MD_COMMONS, ONLY: NOPT, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         LOGICAL, INTENT(IN) :: CENTRET                 
         REAL(KIND=REAL64)

      END SUBROUTINE REMOVE_LINMOM
   ! 1. get random distributed velocities and determine net translation
   DO I=1,NATOMS
      DO J=1,3
         ! random normal-distributed velocities with zero mean and unit standard deviation  
         CALL RAND_NORMAL(1.0D0, 0.0D0, CURRVEL)
         VEL(3*(I-1)+J) = CURRVEL
         ! get contributions to centre of mass and initial momentum
         COM(J) = COM(J) + COORDS(3*(I-1)+J)*MASSES(I)
         P_COM(J) = P_COM(J) + CURRVEL*MASSES(I)
      END DO
   END DO
   COM(1:3) = COM(1:3)/TOTALMASS
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Centre of mass:             ", &
                                COM(1), COM(2), COM(3)
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Momentum of centre of mass: ", &
                                P_COM(1), P_COM(2), P_COM(3)

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