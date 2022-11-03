MODULE MD_SIMULATION
   CONTAINS
      SUBROUTINE ZERO_STEP()
         USE MD_COMMONS, ONLY: NATOMS, MYUNIT, COORDS, EPOT, EKIN, ACC, MASSES
         USE MD_UTILS, ONLY: SET_DERIVED_PARAMS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         IMPLICIT NONE
         INTEGER :: I, J, IDX
         ! set half step and friction params
         CALL SET_DERIVED_PARAMS()
         ! get initial energies
         CALL HIRE_ENERGY_GRAD(3*NATOMS, COORDS, EPOT, ACC)
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ACC(IDX) = -ACC(IDX)/MASSES(I)
            END DO
         END DO
         WRITE(MYUNIT,'(A)') " mdhire> Initial energies - EPOT= ", EPOT, "; EKIN= ", EKIN
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
         USE MD_COMMONS, ONLY: MYUNIT, NDUMPE, HDT, DT, GAMMA, GFRIC, COORDS, VEL, ACC, MASSES, EKIN, EPOT 
         USE RAND_ROUTINES, ONLY: RAND_NORMAL
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         REAL(KIND = REAL64), PARAMETER :: PI = 3.141592653589793D0
         REAL(KIND = REAL64) :: NR1, NR2
         REAL(KIND = REAL64) :: NOISE(NATOMS)
         INTEGER :: I, J, K, IDX
         ! set noise to be used
         DO I=1,NATOMS
            NOISE(I) = DSQRT(TEMP*GAMMA*DT/MASSES(I))
         END DO

         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ! half-step velocity update
               CALL RAND_NORMAL(1.0, 0.0, NR1)
               VEL(IDX) = GFRIC*VEL(IDX) + ACC(IDX)*HDT + NR1*NOISE(I)
               !update full-step coordinates
               COORDS(IDX) = COORDS(IDX) + VEL(IDX)*DT
            END DO
         END DO
         ! get new potential energy and gradient
         CALL HIRE_ENERGY_GRAD(3*NATOMS, COORDS, EPOT, ACC)
         EKIN=0.0D0
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1) + J
               ! update acceleration
               ACC(IDX) = -ACC(IDX)/MASSES(I)
               ! update full step velocity
               CALL RAND_NORMAL(1.0, 0.0, NR2)
               VEL(IDX) = GFRIC*VEL(IDX) + ACC(IDX)*HDT + NR2*NOISE(I)
               ! get kinetic energy
               EKIN = EKIN + MASSES(I)*VEL(IDX)*VEL(IDX)
            END DO
         END DO
         EKIN = 0.5*EKIN
         IF (MOD(CURRSTEP,NDUMPE).EQ.0) THEN
            WRITE(MYUNIT,*) " mdsteps> Completed step ", CURRSTEP
            WRITE(MYUNIT,*) "          Total energy: ", EPOT+EKIN           
            WRITE(MYUNIT,*) "          Kinetic energy: ", EKIN
            WRITE(MYUNIT,*) "          Potential energy: ", EPOT
            WRITE(MYUNIT,*) " -----------------------------------------------"
         END IF
      END SUBROUTINE TAKE_MDSTEP

      SUBROUTINE DUMPDATA(CURRSTEP)
         USE MD_COMMONS, ONLY: NATOMS, XUNIT, EUNIT, DUMPPDBT, NDUMPE, NDUMPP, NDUMPX, COORDS, EKIN, EPOT
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: CURRSTEP
         INTEGER :: I
         CHARACTER(LEN=15) :: JSTRING
         CHARACTER(LEN=30) :: PDBNAME
         !write energies
         WRITE(EUNIT,(I10,3(1X,F15.7))) CURRSTEP, EPOT+EKIN, EPOT, EKIN
         !write coordinate files
         IF (MOD(CURRSTEP,NDUMPX).EQ.0) THEN
            WRITE(XUNIT,*) " Step: ", CURRSTEP
            DO I=1,NATOMS
               WRITE(XUNIT,'(3F15.7)') COORDS(3*I-2), COORDS(3*I-1), COORD(3*I)
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