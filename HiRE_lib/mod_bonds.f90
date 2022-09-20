!> @file
!> Contains MOD_BONDS to handle CG bonds

!> Module containing all routines and variables to calculate the bonding contribution to E and F\n
!> Bonds: E = RK * (R - REQ)**2 with R the current distance
MODULE MOD_BONDS
   USE PREC_HIRE
   !> Number of bonds
   INTEGER :: NBONDS   
   !> Number of unique bond types 
   INTEGER :: NUMBND   
 
   !Bonds: E = RK * (R - REQ)**2 with R the current distance
   !> Force constants for bonds
   REAL(KIND = REAL64), ALLOCATABLE :: RK(:)     
   !> Equilibrium bond lengths
   REAL(KIND = REAL64), ALLOCATABLE :: REQ(:)    
  
   !bonding book-keeping
   !> Bond atom1
   INTEGER, ALLOCATABLE :: IB(:)                 
   !> Bond atom2
   INTEGER, ALLOCATABLE :: JB(:)                 
   !> Type of bond 
   INTEGER, ALLOCATABLE :: ICB(:)                

   CONTAINS
      !> Routine to allocate all required arrays
      SUBROUTINE ALLOC_BONDS()
         CALL DEALLOC_BONDS()
         ALLOCATE(RK(NUMBND), REQ(NUMBND), IB(NBONDS), JB(NBONDS), ICB(NBONDS))
      END SUBROUTINE ALLOC_BONDS
      !> Calculate the energy and force contribution from the bonded terms
      SUBROUTINE ENERGY_BONDS(NOPT, X, F, EBOND)
         USE NAPARAMS, ONLY: SCORE_RNA
         IMPLICIT NONE
      
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EBOND

         REAL(KIND = REAL64) :: XA(3), RIJ, DA, DF
         INTEGER :: JN, I, J, IC
            
         !initialise force and energy
         F(1:NOPT) = 0.0d0
         EBOND = 0.0d0
         DO JN = 1,NBONDS
            I = IB(JN)
            J = JB(JN)
            XA = X((I+1):(I+3))-X((J+1):(J+3)) 
            RIJ = DSQRT(DOT_PRODUCT(XA,XA))    
            IC = ICB(JN)
            DA = RIJ-REQ(IC)
            DF = RK(IC)*DA*SCORE_RNA(1)
            EBOND = EBOND + DF*DA
            DF = (DF+DF)/RIJ
            XA = DF*XA
            F((I+1):(I+3)) = F((I+1):(I+3)) - XA
            F((J+1):(J+3)) = F((J+1):(J+3)) + XA
         END DO

      END SUBROUTINE ENERGY_BONDS

      !> Deallocate all arrays in this module  
      SUBROUTINE DEALLOC_BONDS()
         IF (ALLOCATED(RK)) DEALLOCATE(RK)
         IF (ALLOCATED(REQ)) DEALLOCATE(REQ)
         IF (ALLOCATED(IB)) DEALLOCATE(IB)
         IF (ALLOCATED(JB)) DEALLOCATE(JB)
         IF (ALLOCATED(ICB)) DEALLOCATE(ICB)
      END SUBROUTINE DEALLOC_BONDS

END MODULE MOD_BONDS
