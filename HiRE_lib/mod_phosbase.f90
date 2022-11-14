!> @file
!> Contains MOD_PHOSBASE dealing with phosphate base interactions

!> Module containing all routines and variables to calculate the phosphate-base interaction contributions to E and F\n
!> 
MODULE MOD_PHOSBASE
   USE PREC_HIRE
   !> Put module variables here              

   CONTAINS


      !> Calculate the energy and force contribution from the phosphate-base interactions
      SUBROUTINE ENERGY_PHOSBASE(NOPT, X, F, EPB)
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from phosphate base interaction
         REAL(KIND = REAL64), INTENT(OUT) :: EPB

      END SUBROUTINE ENERGY_PHOSBASE


END MODULE MOD_PHOSBASE
