SUBROUTINE INITIALISE_VEL()
   USE NUMKIND
   USE MD_COMMONS, ONLY: NATOMS, TEMP, COORDS, VEL, MASSES, MYUNIT, EKIN
   USE RAND_ROUTINES, ONLY: RAND_NORMAL
   USE UTILS_VEC, ONLY: CROSSP
   IMPLICIT NONE
   REAL(KIND = REAL64) :: COM(3)      ! centre of mass
   REAL(KIND = REAL64) :: P_COM(3)    ! linear momentum of the centre of mass
   REAL(KIND = REAL64) :: TOTALMASS   ! total mass
   REAL(KIND = REAL64) :: CURRVEL     ! place holder for current velocity
   REAL(KIND = REAL64) :: ANG_MOM(3)  ! angular momentum
   REAL(KIND = REAL64) :: MOI(3,3)    ! moment of inertia
   REAL(KIND = REAL64) :: R(3)        ! position vector
   REAL(KIND = REAL64) :: R2          ! squared length of position vector
   REAL(KIND = REAL64) :: P(3)        ! momentum vector 
   REAL(KIND = REAL64) :: W(3)        ! angular velocity    
   REAL(KIND = REAL64) :: X(3*NATOMS) ! new coordinates, moved to CoM
   REAL(KIND = REAL64) :: WORK(3)     ! workspace array for inversion in LAPACK
   REAL(KIND = REAL64) :: ATX(3)      ! coordinates for current atom
   REAL(KIND = REAL64) :: SCALE       ! scaling to get correct T
   INTEGER :: I, J, K, IDX, PIVOT(3), STAT


   ! get total mass of system
   TOTALMASS = SUM(MASSES)
   ! initialise COM and P_COM
   COM(1:3) = 0.0D0
   P_COM(1:3) = 0.0D0
   ! initialise angular momentum and moment of inertia
   ANG_MOM(1:3) = 0.0D0
   MOI(1:3,1:3) = 0.0D0
   ! initialise angular velocity and kinetic energy
   W(1:3) = 0.0D0
   EKIN = 0.0D0

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
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Angular momentum:           ", ANG_MOM(1:3)
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Angular velocity:           ", W(1:3)

   ! 3. remove rotation and compute kinetic energy
   ! reset angular momentum and linear momentum
   ANG_MOM(1:3) = 0.0D0
   P_COM(1:3) = 0.0D0
   DO I=1,NATOMS
      ATX(1:3) = X((3*I-2):(3*I))
      R = CROSSP(W,ATX)
      DO J=1,3
         IDX = 3*(I-1) + J
         VEL(IDX) = VEL(IDX) - R(J)
         EKIN = EKIN + MASSES(I)*VEL(IDX)**2 ! the factor of 0.5 is omitted here as it is applied later
         P(J) = VEL(IDX)*MASSES(I)
         P_COM(J) = P_COM(J) + P(J)
      END DO
      ANG_MOM = ANG_MOM + CROSSP(ATX,P)
   END DO
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Final momentum CoM:         ", P_COM(1:3)
   WRITE(MYUNIT,'(A,3(F15.6))') " init_velocity> Final angular momentum:     ", ANG_MOM(1:3)

   ! 4. rescale velocities to match temperature
   SCALE = SQRT(DBLE(3*NATOMS)*TEMP/EKIN)
   EKIN = 0.0D0
   DO I=1,NATOMS
      DO J=1,3
        IDX = 3*(I-1) + J
        VEL(IDX) = SCALE*VEL(IDX)
        EKIN = EKIN + MASSES(I)*VEL(IDX)*VEL(IDX)
      END DO
   END DO
   EKIN = 0.5*EKIN
   WRITE(MYUNIT,'(A,F15.6)') " init_velocity> Kinetic energy: ", EKIN
   WRITE(MYUNIT,'(A)') " "
END SUBROUTINE INITIALISE_VEL