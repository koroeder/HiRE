!> @file
!> Contains setup routine for HiRE

!> Setup for HiRE
!> @brief
!>
!> Routine is setting up the HiRE interface\n
!> First the topology and the scale.dat file are read.\n
!> Then the various modules for energy contributions are initialised and the routine checks for restraints.\n
!> Finally the parameter sets for bases and particles are set and filled.
!>
!> @param[in] TOPNAME - topology name
!> @param[in] SCALEDATNAME - scale.dat file name
SUBROUTINE HIRE_SETUP(TOPNAME, SCALEDATNAME)
   USE MOD_INIT, ONLY: TOPOLOGY_INPUT, READ_SCALE_DAT, CREATE_CONSTR, &
                       INIT_FROM_MODS
   IMPLICIT NONE
   CHARACTER(LEN=*), INTENT(IN) :: TOPNAME !Name of topology file 
   CHARACTER(LEN=*), INTENT(IN) :: SCALEDATNAME !Name of scale.dat file 
   
   !1. Call topology reader
   CALL TOPOLOGY_INPUT(TOPNAME)
   
   !2. Read scale.dat file
   CALL READ_SCALE_DAT(SCALEDATNAME)
   
   !3. Initialise various properties from various modules
   CALL INIT_FROM_MODS()
   
   !4. Check whether we have a restraint file, and if so read in its values
   CALL CREATE_CONSTR()

END SUBROUTINE HIRE_SETUP

