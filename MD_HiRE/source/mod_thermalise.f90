MODULE MOD_VEL
   USE NUMKIND
   USE MD_COMMONS, ONLY: NATOMS, MASSES, MYUNIT
   IMPLICIT NONE
   INTEGER :: NTHERMALISE
   INTEGER :: NEQUIL
   REAL(KIND = REAL64), ALLOCATABLE :: TEMPS(:)
   CONTAINS
      
      SUBROUTINE THERMALISE(TINIT, TFINAL, VELT)
         
         ! Get the temperatures to be used in thermalisation
         ! We set the initial temperature to close to zero
         ! We then use equally spaced intervals
         ALLOCATE(TEMPS(NTHERMALISE))
         TEMPS(1) = TINIT
         TEMPS(NTHERMALISE) = TFINAL
         DTEMP = (TFINAL-TINIT)/DBLE(NTHERMALISE-1)
         DO I = 2,NTHERMALISE-1
            TEMPS(I) = TINIT + I*DTEMP
         END DO

         IF (.NOT.VELT) THEN
            IF (TINIT.LT.1.0D-10) THEN
               TINIT = 1.0D-6
            END IF
            CALL INITIALISE_VEL(TINIT)
         END IF

         DO I=1,NTHERMALISE
            TEMP = TEMPS(I)
            CALL THERMALISE_RESCALE_VEL(VEL,TEMP)

            DO J=1,NEQUIL
               ! add if clauses to check whether we want to remove these
               CALL REMOVE_LINMOM(X, VEL, .TRUE.)
               CALL REMOVE_ANGVEL(X, VEL)
               ! Velocity verlet?
               IF (MDMETHOD.EQ.'VV') THEN
                  ! add clause to check when we rescale
                  CALL SCALEVEL(TEMP,VEL)
                  CALL VELOCITY_VERLET(X, VEL, ACC, EPOT)
               ELSE IF (MDMETHOD.EQ.'LD') THEN
                  ! add clause to check when we rescale
                  CALL SCALEVEL_LANGEVIN(TEMP,VEL)
                  CALL LANGEVIN_STEP(TEMP,X, VEL, ACC, EPOT)
               ELSE  
                  WRITE(MYUNIT,*) " thermalise> No valid MD steps detected"
                  STOP                
               END IF
            END DO

         END DO
         
         !from thermalize routine

         ! 1. remove angular momentum and rotations
         ! 2. cycle over thermalisation steps
         !   a) rescale velocity
         !   b) cycle of equilibration steps
         !      i) potentially remove total lin mom and ang mom
         !      ii) potentially rescale velocities
         !      iii) velocity Verlet step
         !   c) report statistics for this

         ! rescaling velocity 2 (a)
         ! 1. get current temperature (kinetic energy divided by #dof)
         ! 2. rescale with scale = sqrt(target temp/current temp)
         !    vel = vel*scale
      END SUBROUTINE THERMALISE


      SUBROUTINE THERMALISE_RESCALE_VEL(VEL,TEMP)
         USE MD_COMMONS, ONLY: NATOMS, MASSES
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(3*NATOMS)
         REAL(KIND=REAL64), INTENT(IN) :: TEMP
         REAL(KIND=REAL64) :: EKIN, CURRTEMP, ATMASS, SCALE
         INTEGER :: DOF, IDX

         EKIN= 0.0D0
         DOF = 3*NATOMS - 6
         DO I=1,NATOMS
            ATMASS = MASSES(I)
            DO J=1,3
               IDX = 3*(I-1)+J
               EKIN = EKIN + VEL(IDX)*VEL(IDX)/ATMASS
            END DO
         END DO
         
         CURRTEMP = EKIN/DBLE(DOF)
         SCALE = DSQRT(TEMP/CURRTEMP)
         VEL(1:3*NATOMS) = VEL(1:3*NATOMS) * SCALE
      END SUBROUTINE THERMALISE_RESCALE_VEL
END MODULE MOD_VEL