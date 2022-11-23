MODULE MOD_THERMALISE
   USE NUMKIND
   USE MD_COMMONS, ONLY: NATOMS, MASSES, MYUNIT
   IMPLICIT NONE
   !> Number of thermalisation steps
   INTEGER :: NTHERMALISE
   !> Number of MD steps to equilibrate per thermalisation step
   INTEGER :: NEQUIL
   !> Frequency of recentering and removing linear momentum, if 0, will not be used
   INTEGER :: NCENTRE
   !> Frequency of removing ang velocity, if 0, will not be used
   INTEGER :: NRMANG
   !> Frequency of rescaling
   INTEGER :: NRESCALE
   !> Switch whether velocities need to be initialised
   LOGICAL :: VELT = .TRUE.
   CONTAINS
      
      SUBROUTINE THERMALISE(TINIT, TFINAL, X, VEL, ACC, EPOT)

         REAL(KIND = REAL64), INTENT(IN) :: TINIT
         REAL(KIND = REAL64), INTENT(IN) :: TFINAL
         REAL(KIND = REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND = REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND = REAL64), INTENT(INOUT) :: ACC(NOPT)
         REAL(KIND = REAL64), INTENT(OUT) :: EPOT                  
         REAL(KIND = REAL64), ALLOCATABLE :: TEMPS(NTHERMALISE)
         REAL(KIND = REAL64) :: DTEMP, EPOT
         ! Get the temperatures to be used in thermalisation
         ! If the initial T is zero, we set it to be very small, but non-zero
         ! We then use equally spaced intervals
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
            VELT = .FALSE.
         END IF

         DO I=1,NTHERMALISE
            TEMP = TEMPS(I)
            CALL THERMALISE_RESCALE_VEL(VEL,TEMP)

            DO J=1,NEQUIL
               IF (NCENTRE.GT.0) THEN
                  IF (MOD(J,NCENTRE).EQ.0) THEN
                     CALL REMOVE_LINMOM(X, VEL, .TRUE.)
                  END IF
               END IF
               IF (NRMANG.GT.0) THEN
                  IF (MOD(J,NRMANG).EQ.0) THEN
                     CALL REMOVE_ANGVEL(X, VEL)
                  END IF
               END IF
               ! Velocity verlet?
               IF (MDMETHOD.EQ.'VV') THEN
                  IF (MOD(J,NRESCALE).EQ.0) THEN
                     CALL SCALEVEL(TEMP,VEL)
                  END IF
                  CALL VELOCITY_VERLET(X, VEL, ACC, EPOT)
               ! Langevin?
               ELSE IF (MDMETHOD.EQ.'LD') THEN
                  IF (MOD(J,NRESCALE).EQ.0) THEN                  
                     CALL SCALEVEL_LANGEVIN(TEMP,VEL)
                  END IF
                  CALL LANGEVIN_STEP(TEMP,X, VEL, ACC, EPOT)
               ELSE  
                  WRITE(MYUNIT,*) " thermalise> No valid MD steps detected"
                  STOP                
               END IF
            END DO

         END DO
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
END MODULE MOD_THERMALISE