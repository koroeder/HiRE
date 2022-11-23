MODULE MD_CALCS
   USE NUMKIND
   CONTAINS
      ! function to get kinetic energy
      REAL(KIND=REAL64) FUNCTION E_KINETIC(VEL) RESULT(EKIN)
         USE MD_COMMONS, ONLY: NOPT, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT)
         INTEGER :: I
         EKIN = 0.0D0
         DO I=1,NOPT
            EKIN = EKIN + VEL(I)*VEL(I)/MASSES(I)
         END DO
         EKIN = 0.5*EKIN
      END FUNCTION E_KINETIC

      ! subroutine to get accelaration form gradient
      SUBROUTINE GET_ACC(GRAD,ACC)
         USE MD_COMMONS, ONLY: NATOMS, MASSES
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: GRAD(3*NATOMS)
         REAL(KIND=REAL64), INTENT(OUT) :: ACC(3*NATOMS)
         INTEGER :: I, J, IDX

         ACC(1:3*NATOMS) = 0.0D0
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ! the HiRE interface changes the sign of GRAD for global optimisation and minimisers
               ! we need to revert that change here for the acceleration
               ACC(IDX) = -GRAD(IDX)/MASSES(I)
            END DO
         END DO
      END SUBROUTINE GET_ACC

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
         WRITE(MYUNIT,'(A,3(F15.6))') " linmom> Centre of mass:             ", &
                                      COM(1), COM(2), COM(3)
         WRITE(MYUNIT,'(A,3(F15.6))') " linmom> Momentum of centre of mass: ", &
                                      PCOM(1), PCOM(2), PCOM(3)
      END SUBROUTINE DETERMINE_LINMOM

      SUBROUTINE REMOVE_LINMOM(X,VEL,CENTRET)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         LOGICAL, INTENT(IN) :: CENTRET                 
         REAL(KIND=REAL64) :: COM(3), PCOM(3)
         INTEGER :: I, J, IDX

         CALL DETERMINE_LINMOM(X,VEL,COM,PCOM)
         TOTALMASS = SUM(MASSES)/3
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J 
               VEL(IDX) = VEL(IDX) - PCOM(J)/TOTALMASS
               X(IDX) = X(IDX) - COM(J)
            END DO
         END DO
      END SUBROUTINE REMOVE_LINMOM


      SUBROUTINE GET_ANGMOM(X,VEL,ANGMOM,W)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         REAL(KIND = REAL64), INTENT(OUT) :: ANGMOM(3)   ! angular momentum
         REAL(KIND = REAL64), INTENT(OUT) :: W(3)        ! angular velocity 
         REAL(KIND = REAL64) :: MOI(3,3)    ! moment of inertia
         REAL(KIND = REAL64) :: R(3)        ! position vector
         REAL(KIND = REAL64) :: R2          ! squared length of position vector
         REAL(KIND = REAL64) :: P(3)        ! momentum vector 
         INTEGER :: I, J, K, IDX

         ! initialise angular momentum and moment of inertia
         ANG_MOM(1:3) = 0.0D0
         MOI(1:3,1:3) = 0.0D0
         ! initialise angular velocity
         W(1:3) = 0.0D0
         DO I=1,NATOMS
            R(1:3) = 0.0D0
            P(1:3) = 0.0D0
            DO J=1,3
               IDX = 3*(I-1) + J
               R(J) = X(3*(I-1) + J)
               P(J) = VEL(IDX) * MASSES(I)
            END DO
            R2 = DOT_PRODUCT(R,R)
            ANG_MOM = ANG_MOM + CROSSP(R, P)
            DO J=1,3
               MOI(J,J) = MOI(J,J) + (R2 - R(J)*R(J))*MASSES(I)
               DO K=J,3
                  MOI(J,K) = MOI(J,K) - R(J)*R(K)*MASSES(I)
                  MOI(K,J) = MOI(J,K)
               END DO
            END DO                        
         END DO

         ! get LU factorisation for moment of inertia
         CALL DGETRF(3,3,MOI,3,PIVOT,STAT)
         IF (STAT.EQ.0) THEN
            ! get inverse matrix for moment of inertia
            CALL DGETRI(3,MOI,3,PIVOT,WORK,3,STAT)
            IF (STAT.NE.0) THEN
              WRITE(MYUNIT,*) " init_velocity> Error in matrix inversion with DGETRI, error code: ", STAT
            END IF
         ELSE 
            WRITE(MYUNIT,*) " init_velocity> Error in LU factorisation with DGETRF, error code: ", STAT
            STOP
         END IF
         ! compute the angular velocity
         DO J=1,3
            W(J) = DOT_PRODUCT(MOI(J,1:3),ANG_MOM(1:3))
         END DO
         WRITE(MYUNIT,'(A,3(F15.6))') " get_angmom> Angular momentum:           ", ANG_MOM(1:3)
         WRITE(MYUNIT,'(A,3(F15.6))') " get_angmom> Angular velocity:           ", W(1:3)      
      END SUBROUTINE GET_ANGMOM

      SUBROUTINE REMOVE_ANGVEL(X,VEL)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         REAL(KIND = REAL64) :: ANGMOM(3)   ! angular momentum
         REAL(KIND = REAL64) :: ANG_MOM(3)  ! angular momentum         
         REAL(KIND = REAL64) :: W(3)        ! angular velocity 
         REAL(KIND = REAL64) :: R(3)        ! position vector
         REAL(KIND = REAL64) :: ATX(3)      ! atomic coords
         REAL(KIND = REAL64) :: P(3)        ! momentum vector 
         INTEGER :: I, J, IDX

         CALL GET_ANGMOM(X,VEL,ANGMOM,W)

         ANG_MOM(1:3) = 0.0D0
         DO I=1,NATOMS
            ATX(1:3) = X((3*I-2):(3*I))
            R = CROSSP(W,ATX)
            DO J=1,3
               IDX = 3*(I-1) + J
               VEL(IDX) = VEL(IDX) - R(J)
               P(J) = VEL(IDX)*MASSES(I)
            END DO
            ANG_MOM = ANG_MOM + CROSSP(ATX,P)
         END DO
         WRITE(MYUNIT,'(A,3(F15.6))') " rm_angvel> Final angular momentum:     ", ANG_MOM(1:3)
      END SUBROUTINE REMOVE_ANGVEL
END MODULE MD_CALCS