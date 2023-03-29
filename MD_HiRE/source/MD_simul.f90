MODULE MD_SIMULATION
   CONTAINS
      SUBROUTINE ZERO_STEP()
         USE MD_COMMONS, ONLY: NATOMS, MYUNIT, COORDS, EPOT, EKIN, ACC, VEL, TEMP, &
                               NOPT, MDMETHOD, TFINAL, TINIT, THERMINIT, RESTARTSIMT, &
                               RESTARTINPF, RESTARTSTEP, RMSDT, THERMINIT2
         USE MD_UTILS, ONLY: SET_DERIVED_PARAMS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MD_CALCS, ONLY: GET_ACC, E_KINETIC
         USE MOD_THERMALISE, ONLY: THERMALISE, THERMALISE2
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
            IF (THERMINIT.OR.THERMINIT2) THEN
               IF (TFINAL.LT.0.0D0) THEN
                  TFINAL = TEMP
               END IF
               WRITE(MYUNIT,*) " mdhire> Thermalisation from ", TINIT, " to ", TFINAL   
               IF (THERMINIT2) THEN 
                  CALL THERMALISE2(TINIT, TFINAL, COORDS, VEL, ACC, EPOT)
               ELSE IF (THERMINIT) THEN
                  CALL THERMALISE(TINIT, TFINAL, COORDS, VEL, ACC, EPOT)
               END IF
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
         CALL FLUSH(MYUNIT)
      END SUBROUTINE ZERO_STEP

      SUBROUTINE RUN_MD()
         USE MD_COMMONS, ONLY: MDSTEPS, RESTARTSTEP, CONTINUESIMT, REXT, TASKID, NREPLICA, MYUNIT
         USE EXCHANGES, ONLY: CURRENT_ORDER, MYCURRENTID, UPDATE_CURR_ORDER
         USE NUMKIND
         IMPLICIT NONE
         INTEGER :: J
         CHARACTER(LEN=10) :: DATECHAR, TIMECHAR, ZONECHAR
         INTEGER :: VALUES(8), ITIME
         REAL(KIND=REAL64) :: DPRAND, CURRTEMP
 
         IF (REXT) THEN
            WRITE(MYUNIT, '(A)') " mdrun> Initialising for REX simulation"
            IF (.NOT.ALLOCATED(CURRENT_ORDER)) ALLOCATE(CURRENT_ORDER(NREPLICA))
            MYCURRENTID = TASKID + 1
            CALL UPDATE_CURR_ORDER()

            CALL DATE_AND_TIME(DATECHAR,TIMECHAR,ZONECHAR,VALUES)
            ITIME = VALUES(6)*60 + VALUES(7)
            CALL SDPRND(ITIME+TASKID)
         END IF
         CALL FLUSH(MYUNIT)
         IF (CONTINUESIMT) THEN
            IF (RESTARTSTEP.LT.MDSTEPS) THEN
               DO J=RESTARTSTEP, MDSTEPS
                  IF (REXT) CALL EXCHANGEREPS(J)
                  CALL TAKE_MDSTEP(J,CURRTEMP)
                  CALL DUMPDATA(J,CURRTEMP)
               END DO 
            END IF
         ELSE
            DO J=1,MDSTEPS
               CALL TAKE_MDSTEP(J,CURRTEMP)
               CALL DUMPDATA(J,CURRTEMP)
               IF (REXT) CALL EXCHANGEREPS(J)
            END DO
         END IF
          IF (REXT) THEN
            IF (.NOT.ALLOCATED(CURRENT_ORDER)) DEALLOCATE(CURRENT_ORDER)
         END IF
      END SUBROUTINE RUN_MD

      SUBROUTINE EXCHANGEREPS(J)
         USE EXCHANGES, ONLY: SELECT_EXCHANGES
         USE MD_COMMONS, ONLY: NREXSTEPS, MDSTEPS, TASKID, MYUNIT
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: J

         ! avoid exchanges towards the end of the simulation
         IF ((MDSTEPS-J).LT.NREXSTEPS) THEN
            RETURN
         ELSE 
            IF (MOD(J,NREXSTEPS).EQ.0) THEN
               WRITE(MYUNIT,'(A)') " "
               WRITE(MYUNIT,'(A)') " rexmd> Attempting replica exchanges"
               CALL SELECT_EXCHANGES()
            END IF
         END IF
      END SUBROUTINE EXCHANGEREPS

      SUBROUTINE TAKE_MDSTEP(CURRSTEP,CURRTEMP)
         USE NUMKIND
         USE MD_COMMONS, ONLY: MYUNIT, NOPT, NATOMS, NDUMPE, HDT, DT, GAMMA, GFRIC, REXT, &
                               COORDS, VEL, ACC, MASSES, EKIN, EPOT, TEMP, MDMETHOD, TEMPUNIT
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         USE MOD_INTEGRATORS, ONLY: VELOCITY_VERLET, SCALEVEL, LANGEVIN_STEP
         USE MD_CALCS, ONLY: CURRENT_T, REMOVE_LINMOM, REMOVE_ANGVEL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64), INTENT(OUT) :: CURRTEMP
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

         CALL REMOVE_LINMOM(COORDS,VEL,.FALSE.)

         IF (MOD(CURRSTEP,50).EQ.0) THEN
            CALL REMOVE_ANGVEL(COORDS,VEL)
         END IF

         IF (MOD(CURRSTEP,NDUMPE).EQ.0) THEN
            CURRTEMP = CURRENT_T(VEL)
            WRITE(MYUNIT,*) " mdsteps> Completed step ", CURRSTEP
            WRITE(MYUNIT,'(A,F12.4)') "           Total energy:        ", EPOT+EKIN           
            WRITE(MYUNIT,'(A,F12.4)') "           Kinetic energy:      ", EKIN
            WRITE(MYUNIT,'(A,F12.4)') "           Potential energy:    ", EPOT
            WRITE(MYUNIT,'(A,F12.4)') "           Current temperature: ", CURRTEMP
            WRITE(MYUNIT,'(A,F12.4)') " --------------------------------------------------"
            IF (REXT) WRITE(TEMPUNIT,'(I10,F12.4)') CURRSTEP, CURRTEMP, TEMP
            CALL FLUSH(MYUNIT)
         END IF
      END SUBROUTINE TAKE_MDSTEP

      SUBROUTINE DUMPDATA(CURRSTEP, CURRTEMP)
         USE MD_COMMONS, ONLY: NATOMS, XUNIT, EUNIT, DUMPPDBT, NDUMPE, NDUMPP, NDUMPX, NDUMPRST, &
                               COORDS, VEL, EKIN, EPOT, ALIGNCONFT, RMSDT, NDUMPR, RUNIT, &
                               NTASKS, TASKID, ELEMENTS, REXT
         USE HIRE_INTERFACE, ONLY: DUMP_PDB
         USE MOD_RESTART, ONLY: WRITE_RST_FILE
         USE MOD_RMSD, ONLY: GET_RMSD 
         USE NUMKIND
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64), INTENT(IN) :: CURRTEMP
         INTEGER :: I
         REAL(KIND = REAL64) :: DIST, RMSD
         CHARACTER(LEN=15) :: JSTRING
         CHARACTER(LEN=30) :: PDBNAME
         CHARACTER(LEN=8) :: TASKSTR

         !write energies
         IF (MOD(CURRSTEP,NDUMPE).EQ.0) THEN
            WRITE(EUNIT,'(I10,3(1X,F15.7))') CURRSTEP, EPOT+EKIN, EPOT, EKIN
         END IF 
         !write coordinate files
         IF (MOD(CURRSTEP,NDUMPX).EQ.0) THEN
            WRITE(XUNIT,'(I6)') NATOMS 
            WRITE(XUNIT,*) " Step: ", CURRSTEP, "Total energy: ", EPOT+EKIN, "Temperature: ", CURRTEMP
            DO I=1,NATOMS
               WRITE(XUNIT,'(A1,2X,3F15.7)') ELEMENTS(I),COORDS(3*I-2), COORDS(3*I-1), COORDS(3*I)
            END DO
            !WRITE(XUNIT,*) "-----------------------------------------------"
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
               IF (NTASKS.EQ.1) THEN 
                  PDBNAME = "mdx_"//ADJUSTL(TRIM(JSTRING))//".pdb"
               ELSE
                  WRITE(TASKSTR,'(I6)') TASKID
                  PDBNAME = "mdx_"//ADJUSTL(TRIM(JSTRING))//"."//TRIM(ADJUSTL(TASKSTR))//".pdb"
               END IF
               CALL DUMP_PDB(3*NATOMS,COORDS,PDBNAME,.TRUE.)

            END IF
         END IF 
         !write restart file
         IF (MOD(CURRSTEP,NDUMPRST).EQ.0) THEN
            IF (.NOT.REXT) CALL WRITE_RST_FILE(CURRSTEP, COORDS, VEL)
         END IF
      END SUBROUTINE DUMPDATA
END MODULE MD_SIMULATION
