MODULE MD_SIMULATION
   CONTAINS
      SUBROUTINE ZERO_STEP()
         USE MD_COMMONS, ONLY: NATOMS, MYUNIT, COORDS, EPOT, EKIN, ACC, VEL, TEMP, &
                               NOPT, MDMETHOD, TFINAL, TINIT, THERMINIT
         USE MD_UTILS, ONLY: SET_DERIVED_PARAMS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC
         USE MOD_THERMALISE, ONLY: THERMALISE
         USE NUMKIND
         IMPLICIT NONE
         REAL(KIND = REAL64) :: GRAD(NOPT)
         INTEGER :: I, J, IDX
         ! set half step and friction params
         CALL SET_DERIVED_PARAMS()
         ! get initial energies
         CALL HIRE_ENERGY_GRAD(3*NATOMS, COORDS, EPOT, GRAD)
         ! get acceleration
         CALL GET_ACC(GRAD,ACC)
         WRITE(MYUNIT,'(2(A,F12.4))') " mdhire> Initial energies - EPOT= ", EPOT, "; EKIN= ", EKIN
         WRITE(MYUNIT, '(A)') " "
         IF (THERMINIT) THEN
            IF (TFINAL.LT.0.0D0) THEN
               TFINAL = TEMP
            END IF
            WRITE(MYUNIT,*) " mdhire> Thermalisation from ", TINIT, " to ", TFINAL   
            CALL THERMALISE(TINIT, TFINAL, COORDS, VEL, ACC, EPOT)
            IF (TEMP.NE.TFINAL) THEN
               WRITE(MYUNIT,*) " mdhire> WARNING: Final T of thermalisation is not the same as simulation T."
            END IF
         ELSE
            WRITE(MYUNIT,'(A)') " mdhire> Calling velocity initialisation"     
            CALL INITIALISE_VEL(TEMP)
         END IF
      END SUBROUTINE ZERO_STEP

      SUBROUTINE RUN_MD()
         USE MD_COMMONS, ONLY: MDSTEPS
         IMPLICIT NONE
         INTEGER :: J

         DO J=1,MDSTEPS
            CALL TAKE_MDSTEP(J)
            CALL DUMPDATA(J)
         END DO
      END SUBROUTINE RUN_MD

      SUBROUTINE TAKE_MDSTEP(CURRSTEP)
         USE NUMKIND
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, NDUMPE, HDT, DT, GAMMA, GFRIC, &
                               COORDS, VEL, ACC, MASSES, EKIN, EPOT, TEMP, MDMETHOD
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MOD_INTEGRATORS, ONLY: VELOCITY_VERLET, SCALEVEL, LANGEVIN_STEP
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64) :: NR1, NR2
         REAL(KIND = REAL64) :: NOISE(NATOMS)
         INTEGER :: I, J, IDX

         ! Velocity verlet?
         IF (MDMETHOD.EQ.'VV') THEN
            CALL VELOCITY_VERLET(COORDS, VEL, ACC, EPOT, EKIN)
            CALL SCALEVEL(TEMP,VEL)
         ! Langevin?
         ELSE IF (MDMETHOD.EQ.'LD') THEN

            CALL LANGEVIN_STEP(TEMP, COORDS, VEL, ACC, EPOT, EKIN)
            !IF (MOD(J,NRESCALE).EQ.0) THEN                  
            !   CALL SCALEVEL_LANGEVIN(TEMP,VEL)
            !END IF
         ELSE  
            WRITE(MYUNIT,*) " thermalise> No valid MD steps detected"
            STOP                
         END IF

         IF (MOD(CURRSTEP,NDUMPE).EQ.0) THEN
            WRITE(MYUNIT,*) " mdsteps> Completed step ", CURRSTEP
            WRITE(MYUNIT,'(A,F12.4)') "           Total energy:     ", EPOT+EKIN           
            WRITE(MYUNIT,'(A,F12.4)') "           Kinetic energy:   ", EKIN
            WRITE(MYUNIT,'(A,F12.4)') "           Potential energy: ", EPOT
            WRITE(MYUNIT,'(A,F12.4)') " --------------------------------------------------"
         END IF
      END SUBROUTINE TAKE_MDSTEP

      SUBROUTINE DUMPDATA(CURRSTEP)
         USE MD_COMMONS, ONLY: NATOMS, XUNIT, EUNIT, DUMPPDBT, NDUMPE, NDUMPP, NDUMPX, COORDS, EKIN, EPOT
         USE HIRE_INTERFACE, ONLY: DUMP_PDB
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         INTEGER :: I
         CHARACTER(LEN=15) :: JSTRING
         CHARACTER(LEN=30) :: PDBNAME
         !write energies
         IF (MOD(CURRSTEP,NDUMPE).EQ.0) THEN
            WRITE(EUNIT,'(I10,3(1X,F15.7))') CURRSTEP, EPOT+EKIN, EPOT, EKIN
         END IF 
         !write coordinate files
         IF (MOD(CURRSTEP,NDUMPX).EQ.0) THEN
            WRITE(XUNIT,*) " Step: ", CURRSTEP
            DO I=1,NATOMS
               WRITE(XUNIT,'(3F15.7)') COORDS(3*I-2), COORDS(3*I-1), COORDS(3*I)
            END DO
            WRITE(XUNIT,*) "-----------------------------------------------"
         END IF
         !write pdb
         IF (DUMPPDBT) THEN
            IF (MOD(CURRSTEP,NDUMPP).EQ.0) THEN
               WRITE(JSTRING, '(I12.12)') CURRSTEP
               PDBNAME = "mdx_"//ADJUSTL(TRIM(JSTRING))//".pdb"
               CALL DUMP_PDB(3*NATOMS,COORDS,PDBNAME,.TRUE.)
            END IF
         END IF         
      END SUBROUTINE DUMPDATA
END MODULE MD_SIMULATION