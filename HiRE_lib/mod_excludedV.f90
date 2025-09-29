!> @file
!> Contains MOD_EXCLV handling the excluded volume contributions

!> Module containing all routines and variables to calculate the excluded volume contributions
MODULE MOD_EXCLV
   USE PREC_HIRE
   USE NBDEFS
   IMPLICIT NONE
   !> Overall scaling of the contribution
   REAL(KIND = REAL64) :: EXCLV_SCALING
   !> OLD PARAMETERS, no longer in use.
   !> Steepness of function
   REAL(KIND = REAL64) :: EXCL_VOL
   !> Barrier height
   REAL(KIND = REAL64) :: BARRIER
   !> Ratio modifying distance penalty 
   REAL(KIND = REAL64) :: RATIO
  
   CONTAINS
      !> Setting the common variables to the correct values
      SUBROUTINE INIT_EXCLV()
         USE NAPARAMS, ONLY: SCORE_RNA
         EXCLV_SCALING = SCORE_RNA(40)
      END SUBROUTINE INIT_EXCLV

      !> Energy and gradient for simpler potential
      !  X----)              (------X
      !  <----><------------><------>
      !   r1        d           r2
      !  <-------------------------->
      !          sqrt(DA2)
      ! r1+r2 is defined in CT2 array, so D is sqrt(DA2) - CT2
      SUBROUTINE ENERGY_EXV(NOPT, X, F, I, J, TI, TJ, DA2, EEXCL)
         INTEGER, INTENT(IN) :: NOPT                   !Number of degrees of freedom
         INTEGER, INTENT(IN) :: I, J                   !indices of grains
         INTEGER, INTENT(IN) :: TI,TJ                  !type of grains
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force
         REAL(KIND = REAL64), INTENT(OUT) :: EEXCL     !energy contribution
         REAL(KIND = REAL64), INTENT(IN) :: DA2        !input squared distance between particles
         INTEGER, PARAMETER :: VPOWER = 12
         REAL(KIND = REAL64) :: DX(3), CT2, D, DA, DINV, DF

         CT2 = NBCT2(TI,TJ)
         DX(1:3) = X(3*I-2:3*I) - X(3*J-2:3*J)
         DA = DSQRT(DA2)
         D = DA - CT2
         DINV = 1/D
         EEXCL = EXCLV_SCALING*(DINV**VPOWER)

         DF = -VPOWER*EEXCL*DINV/DA

         F(3*I-2:3*I) = F(3*I-2:3*I) + DF*DX(1:3)
         F(3*J-2:3*J) = F(3*J-2:3*J) - DF*DX(1:3)
      END SUBROUTINE ENERGY_EXV
END MODULE MOD_EXCLV



