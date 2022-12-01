MODULE MOD_THERMALISE
   USE NUMKIND
   USE MD_COMMONS, ONLY: NATOMS, MASSES, MYUNIT
   IMPLICIT NONE
   !> Number of thermalisation steps
   INTEGER :: NTHERMALISE = 25
   !> Number of MD steps to equilibrate per thermalisation step
   INTEGER :: NEQUIL = 1000
   !> Frequency of recentering and removing linear momentum, if 0, will not be used
   INTEGER :: NCENTRE = 1
   !> Frequency of removing ang velocity, if 0, will not be used
   INTEGER :: NRMANG = 0
   !> Frequency of rescaling
   INTEGER :: NRESCALE = 1
   !> Freuqncy of writing energies to output file
   INTEGER ::  NEQDUMPE = 10
   !> Switch whether velocities need to be initialised
   LOGICAL :: VELT = .TRUE.
   CONTAINS
      
      SUBROUTINE THERMALISE(TINIT, TFINAL, X, VEL, ACC, EPOT)
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE MD_COMMONS, ONLY: NOPT, MDMETHOD, TEMP, EKIN
         USE MD_CALCS
         USE MOD_INTEGRATORS, ONLY: LANGEVIN_STEP, VELOCITY_VERLET, SCALEVEL, SCALEVEL_LANGEVIN
         REAL(KIND = REAL64), INTENT(IN) :: TINIT
         REAL(KIND = REAL64), INTENT(IN) :: TFINAL
         REAL(KIND = REAL64), INTENT(INOUT) :: X(NOPT)
         REAL(KIND = REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND = REAL64), INTENT(INOUT) :: ACC(NOPT)
         REAL(KIND = REAL64), INTENT(OUT) :: EPOT                  
         REAL(KIND = REAL64) :: TEMPS(NTHERMALISE)
         REAL(KIND = REAL64) :: DTEMP, TSMALL, CURRTEMP
         INTEGER :: I, J, EQUNIT
         ! Get the temperatures to be used in thermalisation
         ! If the initial T is zero, we set it to be very small, but non-zero
         ! We then use equally spaced intervals
         TEMPS(1) = TINIT
         TEMPS(NTHERMALISE) = TFINAL
         DTEMP = (TFINAL-TINIT)/DBLE(NTHERMALISE-1)
         DO I = 2,NTHERMALISE-1
            TEMPS(I) = TINIT + (I-1)*DTEMP
         END DO

         !  DO I=1,NTHERMALISE-1
         !     TEMPS(I) = TINIT+ I*DTEMP
         !  END DO

         IF (VELT) THEN
            IF (TINIT.LT.1.0D-10) THEN
               TSMALL = 1.0D-3
               TEMPS(1) = TSMALL
            END IF
            CALL INITIALISE_VEL(TSMALL)
            !  CALL INITIALISE_VEL(TEMPS(1))
            VELT = .FALSE.
         END IF

         CALL FILE_OPEN("md_ethermalisation.log",EQUNIT,.TRUE.)

         WRITE(MYUNIT,'(A,I8,A,I8,A)') " thermalise> Thermalise simulation in ", NTHERMALISE, " steps with ", NEQUIL, " equilibration steps for each new T"
         WRITE(MYUNIT,*) " "
         DO I=1,NTHERMALISE
            TEMP = TEMPS(I)
            WRITE(MYUNIT,'(A,I8,A,F7.3)') " thermalise> Thermalisation step ", I, " with T=", TEMP
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
                  CALL VELOCITY_VERLET(X, VEL, ACC, EPOT, EKIN)
               ! Langevin?
               ELSE IF (MDMETHOD.EQ.'LD') THEN
                  !IF (MOD(J,NRESCALE).EQ.0) THEN                  
                  !   CALL SCALEVEL_LANGEVIN(TEMP,VEL)
                  !END IF
                  CALL LANGEVIN_STEP(TEMP, X, VEL, ACC, EPOT, EKIN)
               ELSE  
                  WRITE(MYUNIT,'(A)') " thermalise> No valid MD steps detected"
                  STOP                
               END IF
               IF (MOD(J,NEQDUMPE).EQ.0) THEN
                  CURRTEMP = CURRENT_T(VEL)
                  WRITE(MYUNIT,'(A,I8)') " thermalise> Completed equilibration step ", J
                  WRITE(MYUNIT,'(A,F12.4)') "             Total energy:        ", EPOT+EKIN           
                  WRITE(MYUNIT,'(A,F12.4)') "             Kinetic energy:      ", EKIN
                  WRITE(MYUNIT,'(A,F12.4)') "             Potential energy:    ", EPOT
                  WRITE(MYUNIT,'(A,F7.4)') "             Current temperature: ", CURRTEMP
                  WRITE(MYUNIT,'(A,F12.4)') " --------------------------------------------------"
                  WRITE(EQUNIT,'(2I10,2(1X,F7.4),3(1X,F15.7))') I, J, TEMP, CURRTEMP, EPOT+EKIN, EPOT, EKIN 
               END IF
            END DO        
         END DO
         CLOSE(EQUNIT)
         WRITE(MYUNIT,*) " thermalise> Completed thermalisation"
         WRITE(MYUNIT,*) " " 
      END SUBROUTINE THERMALISE


      SUBROUTINE THERMALISE_RESCALE_VEL(VEL,TEMP)
         USE MD_COMMONS, ONLY: NATOMS, NOPT
         USE MD_CALCS, ONLY: CURRENT_T
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: VEL(NOPT)
         REAL(KIND=REAL64), INTENT(IN) :: TEMP
         REAL(KIND=REAL64) :: CURRTEMP, SCALE

         CURRTEMP = CURRENT_T(VEL)
         SCALE = DSQRT(TEMP/CURRTEMP)
         WRITE(MYUNIT,'(3(A,1X,F7.3))') " thermalise_rescale> Target T: ", TEMP, ", current T:", CURRTEMP, ", scaling: ", SCALE
         VEL(1:3*NATOMS) = VEL(1:3*NATOMS) * SCALE
      END SUBROUTINE THERMALISE_RESCALE_VEL
END MODULE MOD_THERMALISE