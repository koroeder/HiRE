!> @file
!> Contains MOD_EXCLV handling the excluded volume contributions

!> Module containing all routines and variables to calculate the excluded volume contributions
MODULE MOD_EXCLV
   USE PREC_HIRE
   USE NBDEFS, ONLY: NBCT2, SOFTNESS
   IMPLICIT NONE
   !> Overall scaling of the contribution
   REAL(KIND = REAL64) :: EXCLV_SCALING
 
   CONTAINS
      !> Setting the common variables to the correct values
      SUBROUTINE INIT_EXCLV()
         USE NAPARAMS, ONLY: SCORE_RNA
         EXCLV_SCALING = SCORE_RNA(5)
      END SUBROUTINE INIT_EXCLV

      !> Energy and gradient for simpler potential
      !  X----)              (------X
      !  <----><------------><------>
      !   r1        d           r2
      !  <-------------------------->
      !          sqrt(DA2)
      ! r1+r2 is defined in CT2 array, so D is sqrt(DA2) - CT2
      SUBROUTINE ENERGY_EXV(NOPT, X, F, I, J, TI, TJ, DA2, DCORR, EEXCL)
         INTEGER, INTENT(IN) :: NOPT                   !Number of degrees of freedom
         INTEGER, INTENT(IN) :: I, J                   !indices of grains
         INTEGER, INTENT(IN) :: TI,TJ                  !type of grains
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(INOUT) :: F(NOPT) !force
         REAL(KIND = REAL64), INTENT(OUT) :: EEXCL     !energy contribution
         REAL(KIND = REAL64), INTENT(IN) :: DA2        !input squared distance between particles
         REAL(KIND = REAL64), INTENT(IN) :: DCORR      !correction for neighbouring nucleotides
         INTEGER, PARAMETER :: VPOWER = 12
         !> Numerical safety floor on the squared distance (Å²), so a (near-)coincident pair
         !> cannot produce a literal division by zero. Far below any physically meaningful
         !> CG bead separation - only matters once a pair has already collapsed.
         REAL(KIND = REAL64), PARAMETER :: DA2FLOOR = 1.0D-2
         !> Numerical safety cap on the contact-distance ratio D = (bead-size sum)/DA.
         !> D=1 is nominal contact, so DMAX=2 already represents a severe clash; capping D
         !> keeps EEXCL/DF finite (and boundedly large, not physically "correct") for any
         !> closer approach instead of letting the D**VPOWER power law diverge.
         !> TODO: this is a numerical regularization, not a fitted physical softness -
         !> tune DMAX/DA2FLOOR against this force field's actual energy scale, and prefer
         !> populating NBDEFS::SOFTNESS with real per-type-pair values once available.
         REAL(KIND = REAL64), PARAMETER :: DMAX = 2.0D0
         REAL(KIND = REAL64) :: S
         REAL(KIND = REAL64) :: DX(3), CT2, D, DA, DF, DA2G

         S = SOFTNESS(TI,TJ)
         CT2 = NBCT2(TI,TJ)
         DX(1:3) = X(3*I-2:3*I) - X(3*J-2:3*J)
         DA2G = MAX(DA2, DA2FLOOR)
         DA = DSQRT(DA2G)

         !D = DA - (CT2 - DCORR)
         !DINV = 1/D
         !WRITE(*,*) I,J,D,DINV, EEXCL
         !EEXCL = EXCLV_SCALING*(DINV**VPOWER)

         D = MIN(((CT2-S)/DA), DMAX)

         EEXCL = EXCLV_SCALING*(DCORR*D**VPOWER)
         DF = VPOWER*EEXCL/DA2G

         F(3*I-2:3*I) = F(3*I-2:3*I) + DF*DX(1:3)
         F(3*J-2:3*J) = F(3*J-2:3*J) - DF*DX(1:3)
      END SUBROUTINE ENERGY_EXV
END MODULE MOD_EXCLV



