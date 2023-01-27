MODULE MD_SIMULATION
   CONTAINS
      SUBROUTINE ZERO_STEP()
         USE MD_COMMONS, ONLY: NATOMS, MYUNIT, COORDS, EPOT, EKIN, ACC, VEL, TEMP, &
                               NOPT, MDMETHOD, TFINAL, TINIT, THERMINIT, RESTARTSIMT, &
                               RESTARTINPF, RESTARTSTEP, RMSDT
         USE MD_UTILS, ONLY: SET_DERIVED_PARAMS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC, E_KINETIC
         USE MOD_THERMALISE, ONLY: THERMALISE
         USE MOD_RESTART, ONLY: READ_RST_FILE
         USE MOD_RMSD, ONLY: SET_REF
         USE NUMKIND
         IMPLICIT NONE
         REAL(KIND = REAL64) :: GRAD(NOPT)
         REAL(KIND = REAL64) :: TNEW
         INTEGER :: I, J, IDX
         ! set half step and friction params
         CALL SET_DERIVED_PARAMS()
         IF (RESTARTSIMT) THEN
            CALL READ_RST_FILE(RESTARTINPF, RESTARTSTEP, TNEW, COORDS, VEL)
            IF (TNEW.NE.TEMP) THEN
               WRITE(MYUNIT,*) " mdhire> WARNING - Temperature in restart file of ", TNEW, " does not match simulation T of ", TEMP
            END IF
            EKIN = E_KINETIC(VEL)
         END IF

         ! set reference for RMSD
         IF (RMSDT) THEN
            CALL SET_REF(NATOMS,COORDS)
         END IF

         ! get initial energies
         CALL HIRE_ENERGY_GRAD(3*NATOMS, COORDS, EPOT, GRAD)
         ! get acceleration
         CALL GET_ACC(GRAD,ACC)

         IF (.NOT.RESTARTSIMT) THEN
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
         END IF
         WRITE(MYUNIT,'(2(A,F12.4))') " mdhire> Initial energies - EPOT= ", EPOT, "; EKIN= ", EKIN
         WRITE(MYUNIT, '(A)') " "
      END SUBROUTINE ZERO_STEP

      SUBROUTINE RUN_MD()
         USE MD_COMMONS, ONLY: MDSTEPS, RESTARTSTEP, CONTINUESIMT, REXT, TASKID, NREPLICA
         USE EXCHANGES, ONLY: CURRENT_ORDER, MYCURRENTID, UPDATE_CURR_ORDER
         IMPLICIT NONE
         INTEGER :: J

         IF (REXT) THEN
            IF (.NOT.ALLOCATED(CURRENT_ORDER)) ALLOCATE(CURRENT_ORDER(NREPLICA))
            MYCURRENTID = TASKID + 1
            CALL UPDATE_CURR_ORDER()
         END IF
         IF (CONTINUESIMT) THEN
            IF (RESTARTSTEP.LT.MDSTEPS) THEN
               DO J=RESTARTSTEP, MDSTEPS
                  IF (REXT) CALL EXCHANGEREPS(J)
                  CALL TAKE_MDSTEP(J)
                  CALL DUMPDATA(J)
               END DO 
            END IF
         ELSE
            DO J=1,MDSTEPS
               IF (REXT) CALL EXCHANGEREPS(J)
               CALL TAKE_MDSTEP(J)
               CALL DUMPDATA(J)
            END DO
         END IF
          IF (REXT) THEN
            IF (.NOT.ALLOCATED(CURRENT_ORDER)) DEALLOCATE(CURRENT_ORDER)
         END IF
      END SUBROUTINE RUN_MD

      SUBROUTINE EXCHANGEREPS(J)
         USE EXCHANGES, ONLY: SELECT_EXCHANGES
         USE MD_COMMONS, ONLY: NREXSTEPS, MDSTEPS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: J

         ! avoid exchanges towards the end of the simulation
         IF ((MDSTEPS-J).LT.NREXSTEPS) THEN
            RETURN
         ELSE 
            IF (MOD(J,NREXSTEPS).EQ.0) THEN
               CALL SELECT_EXCHANGES()
            END IF
         END IF
      END SUBROUTINE EXCHANGEREPS

      SUBROUTINE TAKE_MDSTEP(CURRSTEP)
         USE NUMKIND
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, NDUMPE, HDT, DT, GAMMA, GFRIC, &
                               COORDS, VEL, ACC, MASSES, EKIN, EPOT, TEMP, MDMETHOD
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MOD_INTEGRATORS, ONLY: VELOCITY_VERLET, SCALEVEL, LANGEVIN_STEP
         USE MD_CALCS, ONLY: CURRENT_T
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64) :: NR1, NR2
         REAL(KIND = REAL64) :: CURRTEMP
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
            CURRTEMP = CURRENT_T(VEL)
            WRITE(MYUNIT,*) " mdsteps> Completed step ", CURRSTEP
            WRITE(MYUNIT,'(A,F12.4)') "           Total energy:        ", EPOT+EKIN           
            WRITE(MYUNIT,'(A,F12.4)') "           Kinetic energy:      ", EKIN
            WRITE(MYUNIT,'(A,F12.4)') "           Potential energy:    ", EPOT
            WRITE(MYUNIT,'(A,F12.4)') "           Current temperature: ", CURRTEMP
            WRITE(MYUNIT,'(A,F12.4)') " --------------------------------------------------"
         END IF
      END SUBROUTINE TAKE_MDSTEP

      SUBROUTINE DUMPDATA(CURRSTEP)
         USE MD_COMMONS, ONLY: NATOMS, XUNIT, EUNIT, DUMPPDBT, NDUMPE, NDUMPP, NDUMPX, NDUMPRST, &
                               COORDS, VEL, EKIN, EPOT, ALIGNCONFT, RMSDT, NDUMPR, RUNIT
         USE HIRE_INTERFACE, ONLY: DUMP_PDB
         USE MOD_RESTART, ONLY: WRITE_RST_FILE
         USE MOD_RMSD, ONLY: GET_RMSD 
         USE NUMKIND
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         INTEGER :: I
         REAL(KIND = REAL64) :: DIST, RMSD
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
         ! get RMSD
         IF (RMSDT) THEN
            IF (MOD(CURRSTEP,NDUMPR).EQ.0) THEN
               CALL GET_RMSD(NATOMS, COORDS, DIST, RMSD, ALIGNCONFT)
               WRITE(RUNIT,'(I10,2(1X,F12.4))') CURRSTEP, DIST, RMSD
            END IF
         END IF
         !write pdb
         IF (DUMPPDBT) THEN
            IF (MOD(CURRSTEP,NDUMPP).EQ.0) THEN
               WRITE(JSTRING, '(I12.12)') CURRSTEP
               PDBNAME = "mdx_"//ADJUSTL(TRIM(JSTRING))//".pdb"
               CALL DUMP_PDB(3*NATOMS,COORDS,PDBNAME,.TRUE.)
            END IF
         END IF 
         !write restart file
         IF (MOD(CURRSTEP,NDUMPRST).EQ.0) THEN
            CALL WRITE_RST_FILE(CURRSTEP, COORDS, VEL)
         END IF
      END SUBROUTINE DUMPDATA
END MODULE MD_SIMULATION