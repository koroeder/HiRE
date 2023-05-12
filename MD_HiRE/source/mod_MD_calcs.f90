MODULE MD_CALCS
   USE NUMKIND
   USE UTILS_VEC
   CONTAINS
      ! function to get kinetic energy
      REAL(KIND=REAL64) FUNCTION E_KINETIC(VEL) RESULT(EKIN)
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT)
         INTEGER :: I, J, IDX
         EKIN = 0.0D0
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               EKIN = EKIN + VEL(IDX)*VEL(IDX)*MASSES(I)
            END DO
         END DO
         EKIN = 0.5*EKIN
      END FUNCTION E_KINETIC

      REAL(KIND=REAL64) FUNCTION CURRENT_T(VEL) RESULT(CURRTEMP)
         USE MD_COMMONS, ONLY: NOPT
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT)
         REAL(KIND=REAL64) :: EKIN
         EKIN = E_KINETIC(VEL)
         CURRTEMP = 2.0D0*EKIN/DBLE(NOPT)         
      END FUNCTION CURRENT_T

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

      ! subroutine to get centre of mass and its momentum
      SUBROUTINE DETERMINE_LINMOM(X,VEL,COM,PCOM)
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, MASSES
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         REAL(KIND=REAL64), INTENT(OUT) :: COM(3)
         REAL(KIND=REAL64), INTENT(OUT) :: PCOM(3)
         INTEGER :: I, J, IDX
         REAL(KIND=REAL64) :: TOTALMASS

         COM(1:3) = 0.0D0
         PCOM(1:3) = 0.0D0
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               COM(J) = COM(J) + X(IDX)*MASSES(I)
               PCOM(J) = PCOM(J) + VEL(IDX)*MASSES(I)              
            END DO
         END DO
         TOTALMASS = SUM(MASSES)
         COM = COM/TOTALMASS
         ! WRITE(MYUNIT,'(A,3(F15.6))') " linmom> Centre of mass:             ", &
         !                             COM(1), COM(2), COM(3)
         ! WRITE(MYUNIT,'(A,3(F15.6))') " linmom> Momentum of centre of mass: ", &
         !                             PCOM(1), PCOM(2), PCOM(3)
      END SUBROUTINE DETERMINE_LINMOM

      SUBROUTINE REMOVE_LINMOM(X,VEL,CENTRET)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         LOGICAL, INTENT(IN) :: CENTRET                 
         REAL(KIND=REAL64) :: COM(3), PCOM(3)
         REAL(KIND = REAL64) :: TOTALMASS
         INTEGER :: I, J, IDX

         CALL DETERMINE_LINMOM(X,VEL,COM,PCOM)
         TOTALMASS = SUM(MASSES)
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J 
               VEL(IDX) = VEL(IDX) - PCOM(J)/TOTALMASS
               X(IDX) = X(IDX) - COM(J)
            END DO
         END DO
      END SUBROUTINE REMOVE_LINMOM


      SUBROUTINE GET_ANGMOM(X,VEL,ANGMOM,W)
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, MASSES
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: VEL(NOPT) 
         REAL(KIND = REAL64), INTENT(OUT) :: ANGMOM(3)   ! angular momentum
         REAL(KIND = REAL64), INTENT(OUT) :: W(3)        ! angular velocity 
         REAL(KIND = REAL64) :: MOI(3,3)    ! moment of inertia
         REAL(KIND = REAL64) :: R(3)        ! position vector
         REAL(KIND = REAL64) :: R2          ! squared length of position vector
         REAL(KIND = REAL64) :: P(3)        ! momentum vector 
         REAL(KIND = REAL64) :: TOTALMASS
         REAL(KIND = REAL64) :: WORK(3)
         INTEGER :: I, J, K, IDX, PIVOT(3), STAT

         ! initialise angular momentum and moment of inertia
         ANGMOM(1:3) = 0.0D0
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
            ANGMOM = ANGMOM + CROSSP(R, P)
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
            W(J) = DOT_PRODUCT(MOI(J,1:3),ANGMOM(1:3))
         END DO
         ! WRITE(MYUNIT,'(A,3(F15.6))') " get_angmom> Angular momentum:           ", ANGMOM(1:3)
         ! WRITE(MYUNIT,'(A,3(F15.6))') " get_angmom> Angular velocity:           ", W(1:3)      
      END SUBROUTINE GET_ANGMOM


      SUBROUTINE GET_ANGMOM2(X,VEL,AMOM_COM,W)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES, TOTALMASS, TMASSINV, MYUNIT
         USE UTILS_VEC, ONLY: CROSSP
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT)          
         REAL(KIND=REAL64) :: COM(3)        ! centre of mass
         REAL(KIND=REAL64) :: PCOM(3)       ! linear momentum of com
         REAL(KIND=REAL64) :: VCOM(3)       ! velocity of com        
         REAL(KIND=REAL64) :: AMOM_COM(3)   ! ang momentum of com
         REAL(KIND=REAL64) :: MOI(3,3)      ! moment of inertia tensor
         REAL(KIND=REAL64) :: W(3)          ! angular velocity       
         REAL(KIND=REAL64) :: XX, XY, XZ, YY, YZ, ZZ, XYZ(3), ATMASS
         REAL(KIND = REAL64) :: WORK(3)
         INTEGER :: I, J, IDX, PIVOT(3), STAT

         CALL DETERMINE_LINMOM(X,VEL,COM,PCOM)

         VCOM(1:3) = PCOM(1:3)*TMASSINV

         AMOM_COM(1:3) = 0.0D0

         DO I=1,NATOMS
            IDX = 3*I-2
            AMOM_COM(1:3) = AMOM_COM(1:3) + CROSSP(X(IDX:IDX+2),VEL(IDX:IDX+2))*MASSES(I)
         END DO

         AMOM_COM(1:3) = AMOM_COM(1:3) - CROSSP(COM,VCOM)*TOTALMASS

         XX = 0.0D0
         XY = 0.0D0
         XZ = 0.0D0
         YY = 0.0D0
         YZ = 0.0D0
         ZZ = 0.0D0

         DO I=1,NATOMS
            XYZ(1) = X(3*I-2) + COM(1)
            XYZ(2) = X(3*I-1) + COM(2)
            XYZ(3) = X(3*I) + COM(3)  
            ATMASS = MASSES(I)
            XX = XX + XYZ(1) * XYZ(1) * ATMASS
            XY = XY + XYZ(1) * XYZ(2) * ATMASS
            XZ = XZ + XYZ(1) * XYZ(3) * ATMASS
            YY = YY + XYZ(2) * XYZ(2) * ATMASS
            YZ = YZ + XYZ(2) * XYZ(3) * ATMASS
            ZZ = ZZ + XYZ(3) * XYZ(3) * ATMASS
         END DO

         MOI(1,1) = YY + ZZ 
         MOI(2,2) = XX + ZZ
         MOI(3,3) = XX + YY
         MOI(1,2) = -XY
         MOI(2,1) = -XY
         MOI(1,3) = -XZ
         MOI(3,1) = -XZ
         MOI(2,3) = -YZ
         MOI(3,2) = -YZ

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
            W(J) = DOT_PRODUCT(MOI(J,1:3),AMOM_COM(1:3))
         END DO

      END SUBROUTINE GET_ANGMOM2


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
         ! WRITE(MYUNIT,'(A,3(F15.6))') " rm_angvel> Final angular momentum:     ", ANG_MOM(1:3)
      END SUBROUTINE REMOVE_ANGVEL

      SUBROUTINE REMOVE_COM_MOTIONS(X, VEL)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, MASSES, TMASSINV
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         REAL(KIND=REAL64) :: COM(3)        ! centre of mass
         REAL(KIND=REAL64) :: PCOM(3)       ! linear momentum of com
         REAL(KIND=REAL64) :: VCOM(3)       ! linear momentum of com         
         REAL(KIND=REAL64) :: ANGMOM(3)     ! angular momentum         
         REAL(KIND=REAL64) :: W(3)          ! angular velocity
         REAL(KIND=REAL64) :: XYZ(3)
         INTEGER :: I, J, IDX

         CALL DETERMINE_LINMOM(X,VEL,COM,PCOM)
         CALL GET_ANGMOM2(X,VEL,ANGMOM,W)

         ! get velocity of CoM
         VCOM(1:3) = PCOM(1:3)*TMASSINV

         ! remove translation
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               VEL(IDX) = VEL(IDX) - VCOM(J)
            END DO
         END DO

         ! stop rotation of CoM
         DO I=1,NATOMS
            DO J=1,3
               XYZ(J) = X(3*(I-1)+J) - COM(J)
            END DO
            VEL(3*I-2) = VEL(3*I-2) - W(2)*XYZ(3) + W(3)*XYZ(2)
            VEL(3*I-1) = VEL(3*I-1) - W(3)*XYZ(1) + W(1)*XYZ(3)
            VEL(3*I)   = VEL(3*I)   - W(1)*XYZ(2) + W(2)*XYZ(1)
         END DO
         CALL GET_ANGMOM(X,VEL,ANGMOM,W)
         
      END SUBROUTINE REMOVE_COM_MOTIONS

      SUBROUTINE CHECKDISTANCE(X,VEL)
         USE MD_COMMONS, ONLY: NATOMS, DIST2REF
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(3*NATOMS)  
         INTEGER :: J
         REAL(KIND=REAL64) :: DIST2, XYZ(3)

         DO J=1,NATOMS
            XYZ = X(3*J-2:3*J)
            DIST2 = DSQ(XYZ(1),XYZ(2),XYZ(3))
            IF (OUTSIDESPHERE(DIST2,DIST2REF)) THEN
               CALL REFLECT_PARTICLE(J, X, VEL)
            END IF
         END DO
      END SUBROUTINE CHECKDISTANCE

      SUBROUTINE REFLECT_PARTICLE(IDMIN, X, VEL)
         USE MD_COMMONS, ONLY: NATOMS
         USE UTILS_VEC
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: IDMIN
         REAL(KIND=REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(3*NATOMS)
         REAL(KIND=REAL64) :: COSB, CURRVEL(3), XYZ(3)

         XYZ(1:3) = VEC_NORMED(X(3*IDMIN-2:3*IDMIN))
         CURRVEL(1:3) = VEC_NORMED(VEL(3*IDMIN-2:3*IDMIN))
         COSB = DOT_PRODUCT(XYZ,CURRVEL)
         CURRVEL(1:3) = CURRVEL(1:3) * (1.0 - 2*COSB)
         VEL(3*IDMIN-2:3*IDMIN) = CURRVEL(1:3)

      END SUBROUTINE REFLECT_PARTICLE


      PURE REAL(KIND=REAL64) FUNCTION DSQ(X, Y, Z)
         REAL(KIND=REAL64), INTENT(IN) :: X, Y, Z
         DSQ = X**2 + Y**2 + Z**2
      END FUNCTION DSQ

      LOGICAL FUNCTION OUTSIDESPHERE(DSQ,DIST2REF) 
         REAL(KIND=REAL64), INTENT(IN) :: DSQ, DIST2REF
         IF (DSQ.LT.DIST2REF) THEN
            OUTSIDESPHERE = .FALSE.
         ELSE 
            OUTSIDESPHERE = .TRUE.
         END IF
      END FUNCTION OUTSIDESPHERE
END MODULE MD_CALCS