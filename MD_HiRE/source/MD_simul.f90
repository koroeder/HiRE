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

END MODULE MD_SIMULATION