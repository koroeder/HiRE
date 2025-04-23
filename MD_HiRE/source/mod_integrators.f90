MODULE MOD_INTEGRATORS
   USE NUMKIND

   CONTAINS

      !> Velocity-verlet algorithm
      !> @brief
      !>
      !> Velocity Verlet algorithm\n
      !> 1. \f( \mathbf{v}(t+\frac{1}{2}\Delta t)\,=\,\mathbf{v}(t)\,+\,\frac{1}{2}\mathbf{a}(t)\Delta t\f) \n
      !> 2. \f( \mathbf{x}(t+\Delta t)\,=\,\mathbf{x}(t)\,+\,\mathbf{v}(t+\frac{1}{2}\Delta t)\Delta t\f) \n
      !> 3. \f( \mathbf{a}(t+\Delta t)\f) from potential\n
      !> 4. \f( \mathbf{V}(t+\Delta t)\,=\,\mathbf{v}(t+\frac{1}{2}\Delta t)\,+\,\frac{1}{2}\mathbf{a}(t+\Delta t)\Delta t\f) \n      
   
      SUBROUTINE VELOCITY_VERLET(X, VEL, ACC, EPOT, EKIN)
         USE MD_COMMONS, ONLY: NATOMS, NOPT, HDT, DT, DETECTHIGHFORCES, WRITINGHIGHFORCES
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC, E_KINETIC
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: ACC(NOPT)
         REAL(KIND=REAL64), INTENT(OUT) :: EPOT, EKIN 
         REAL(KIND=REAL64) :: GRAD(NOPT)

         ! calculate velocity at half step
         VEL(1:NOPT) = VEL(1:NOPT) + ACC(1:NOPT)*HDT
         ! calculate new coordinates
         X(1:NOPT) = X(1:NOPT) + VEL(1:NOPT)*DT
         ! calculate new gradient
         CALL HIRE_ENERGY_GRAD(NOPT, X, EPOT, GRAD, .FALSE.)
         !debugging function for writing trajectories with high forces
         IF (DETECTHIGHFORCES.AND..NOT.(WRITINGHIGHFORCES)) THEN
            CALL HIGH_FORCE_TEST(3*NATOMS,GRAD)          
         END IF
         ! get acceleration form gradient
         CALL GET_ACC(GRAD,ACC)
         ! calculate full step velocity
         VEL(1:NOPT) = VEL(1:NOPT) + ACC(1:NOPT)*HDT
         ! get kinetic energy
         EKIN = E_KINETIC(VEL)
      END SUBROUTINE VELOCITY_VERLET

      !> Velocity-verlet algorithm assuming the acceleration is only position dependent
      SUBROUTINE VELOCITY_VERLET2(X, VEL, ACC, EPOT, EKIN)
         USE MD_COMMONS, ONLY: NATOMS, NOPT, HDT, DT
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC, E_KINETIC
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: ACC(NOPT)
         REAL(KIND=REAL64), INTENT(OUT) :: EPOT, EKIN 
         REAL(KIND=REAL64) :: GRAD(NOPT), OLDACC(NOPT)
 
         ! calculate new coordinates
         X(1:NOPT) = X(1:NOPT) + VEL(1:NOPT)*DT + 0.5*ACC(1:NOPT)*DT*DT
         ! calculate new gradient
         CALL HIRE_ENERGY_GRAD(NOPT, X, EPOT, GRAD, .FALSE.)
         ! save old acceleration
         OLDACC(1:NOPT) = ACC(1:NOPT)
         ! get new acceleration form gradient
         CALL GET_ACC(GRAD,ACC)
         ! calculate full step velocity
         VEL(1:NOPT) = VEL(1:NOPT) + (OLDACC(1:NOPT) + ACC(1:NOPT))*HDT
          ! get kinetic energy
         EKIN = E_KINETIC(VEL)
      END SUBROUTINE VELOCITY_VERLET2      

      SUBROUTINE LANGEVIN_STEP(TEMP,X,VEL,ACC,EPOT, EKIN)
         USE MD_COMMONS, ONLY: NATOMS, NOPT, HDT, DT, GFRIC, GAMMA, MASSES, &
                              DETECTHIGHFORCES, WRITINGHIGHFORCES
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC, E_KINETIC
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: TEMP      
         REAL(KIND=REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND=REAL64), INTENT(INOUT) :: ACC(NOPT)
         REAL(KIND=REAL64), INTENT(OUT) :: EPOT, EKIN        
         REAL(KIND=REAL64) :: GRAD(NOPT)
         REAL(KIND = REAL64) :: NR1, NR2
         REAL(KIND = REAL64) :: NOISE(NATOMS)
         INTEGER :: I, J, IDX 

         ! for random number generation: normal distribution with mean=0, std=1
         ! see Allen and Tildesley, "Computer Simulation of Liquids" (2nd ed.), p. 384 

         DO I=1,NATOMS
            NOISE(I) = DSQRT(TEMP*GAMMA*DT/MASSES(I))
         END DO

         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ! half-step velocity update
               CALL RAND_NORMAL(1.0D0, 0.0D0, NR1)
               VEL(IDX) = GFRIC*VEL(IDX) + ACC(IDX)*HDT + NR1*NOISE(I)
               !update full-step coordinates
               X(IDX) = X(IDX) + VEL(IDX)*DT
            END DO
         END DO
         ! get new potential energy and gradient
         CALL HIRE_ENERGY_GRAD(3*NATOMS, X, EPOT, GRAD, .FALSE.)        
         !debugging function for writing trajectories with high forces
         IF (DETECTHIGHFORCES.AND..NOT.(WRITINGHIGHFORCES)) THEN
            CALL HIGH_FORCE_TEST(3*NATOMS,GRAD)          
         END IF
         ! get acceleration
         CALL GET_ACC(GRAD,ACC)
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ! update full step velocity
               CALL RAND_NORMAL(1.0D0, 0.0D0, NR2)
               VEL(IDX) = GFRIC*VEL(IDX) + ACC(IDX)*HDT + NR2*NOISE(I)
            END DO
         END DO
         ! get kinetic energy
         EKIN = E_KINETIC(VEL)
      END SUBROUTINE LANGEVIN_STEP

      SUBROUTINE SCALEVEL(TEMP,VEL)
         USE MD_COMMONS, ONLY: NOPT
         USE MD_CALCS, ONLY: E_KINETIC
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: TEMP         
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         REAL(KIND=REAL64) :: EKIN, CURRTEMP, SCALE
      
         EKIN = E_KINETIC(VEL)
         ! dof for velocity are 3N, see Allen and Tildesley, "Computer Simulation of Liquids" (2nd ed.), p. 131 
         CURRTEMP = 2.0D0*EKIN/NOPT
         SCALE = DSQRT(TEMP/CURRTEMP)
         VEL(1:NOPT) = VEL(1:NOPT)*SCALE
      END SUBROUTINE SCALEVEL

      SUBROUTINE SCALEVEL_LANGEVIN(TEMP,VEL)
         USE MD_COMMONS, ONLY: NOPT, NATOMS, LANGEVINSCALE, MASSES
         USE MD_CALCS, ONLY: E_KINETIC
         USE RAND_ROUTINES, ONLY: RAND_NORMAL    
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: TEMP         
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT) 
         REAL(KIND=REAL64) :: SCALE, NR
         INTEGER :: I, J, IDX

         SCALE = (1.0D0 - LANGEVINSCALE)*(1.0D0 + LANGEVINSCALE)*TEMP
         DO I=1,NATOMS
            DO J=1,3 
               IDX = 3*(I -1) + J
               CALL RAND_NORMAL(1.0D0, 0.0D0, NR)
               SCALE = SCALE/MASSES(I)*NR
               VEL(IDX) = VEL(IDX)* LANGEVINSCALE + SQRT(SCALE)
            END DO
         END DO  
      END SUBROUTINE SCALEVEL_LANGEVIN

      SUBROUTINE HIGH_FORCE_TEST(NOPT,GRAD)
         USE MD_COMMONS, ONLY: WRITINGHIGHFORCES, FORCETHRESHOLD, HIGHFUNIT
         USE FILE_UTILS, ONLY: FILE_OPEN
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND = REAL64), INTENT(IN) :: GRAD(NOPT)
         REAL(KIND = REAL64) :: RMS

         RMS=MAX(SQRT(SUM(GRAD(1:NOPT)**2)/(NOPT)), 1.0D-100 )

         IF (RMS.GT.FORCETHRESHOLD) THEN
            WRITINGHIGHFORCES = .TRUE.
            WRITE(*,*) " *** High force detected *** "
            WRITE(*,'(2(A,E10.5))') " RMS of ", RMS, " is greater than force threshold of ", FORCETHRESHOLD
            CALL FILE_OPEN("traj_highforces.xyz",HIGHFUNIT,.FALSE.)
         END IF

      END SUBROUTINE HIGH_FORCE_TEST

END MODULE MOD_INTEGRATORS