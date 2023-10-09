!> @file
!> Contains NBDEFS module that handles the variables for the non-bonded terms

!> The module contains the definitions and parameter values for the non-bonded terms.\n
!> It also contains the setup to populate the non-bonded scoring matrix and a function to obtain values from it.
MODULE NBDEFS
   USE PREC_HIRE
   ! Non-bonded variables
   ! NTYPES is the number of particle types, found in the .top file
   ! While the number can potentially change, we have 15 types implemented 
   ! Currently, it's taken to be :
   ! 1: C5*  2: O5*  3: P  4: CA  5: CY
   ! 6,7: G1,2  8,9: A1,2  10: U1  11: C1
   ! 12: D 13: MG 14: NA  15:CL
   ! Update: we have now 16 types, there is CY for RNA and for DNA
   ! 1: C5*  2: O5*  3: P  4: CA  5: CY (RNA), 6: CY (DNA),
   ! 7,8: G1,2  9,10: A1,2  11: U1  12: C1
   ! 13: D 14: MG 15: NA  16: CL
   !> Number of different particles currently implemented
   INTEGER, PARAMETER  :: NTYPES = 16              
   !> Size paramters for general beads
   REAL(KIND = REAL64) :: NBCT2GEN   = 3.6D0
   !> Size paramters for C4 beads   
   REAL(KIND = REAL64) :: NBCT2C4    = 4.0D0
   !> Size paramters for C4 beads 
   REAL(KIND = REAL64) :: NBCT2CY    = 4.0D0 
   !> Size paramters for CY beads 
   REAL(KIND = REAL64) :: NBCT2PP    = 4.0D0
   !> Size paramters for PP beads 
   REAL(KIND = REAL64) :: NBCT2BASE  = 3.2D0
   !> Size paramters for DUMMY beads 
   REAL(KIND = REAL64) :: NBCT2DUMMY = 8.0D0
   !> Effective radii for general beads
   REAL(KIND = REAL64) :: NBRADGEN = 4.0D0
   !> Effective radii for large beads
   REAL(KIND = REAL64) :: NBRADLARGE = 5.0D0
   
   !> Array with 1-4 coefficients
   REAL(KIND = REAL64) :: NBCOEF(NTYPES,NTYPES)
   !> Array with bead sizes
   REAL(KIND = REAL64) :: NBCT2(NTYPES,NTYPES)
   !> Array with effective radii
   REAL(KIND = REAL64) :: NBSCORE(NTYPES,NTYPES)
   !> Charges of particle types, these charges should never change!
   REAL(KIND = REAL64), DIMENSION(NTYPES), PARAMETER ::  &
                         CHRG = (/0.0, 0.0,-1.0, 0.0, 0.0, 0.0, &
                                  0.0, 0.0, 0.0, 0.0, 0.0, &
                                  0.0, 0.0, 2.0, 1.0,-1.0 /) 
                                               
   CONTAINS
      !> Routine to populate all NB arrays
      SUBROUTINE SET_NBPARAMS()
         !generic parameters
         NBCOEF(1:NTYPES,1:NTYPES) = 1.0D0
         NBCT2(1:NTYPES,1:NTYPES) = NBCT2GEN
         NBSCORE(1:NTYPES,1:NTYPES) = NBRADGEN
         !size adjustments
         ! C4 are bigger
         NBCT2(4,1:NTYPES) = NBCT2C4
         NBCT2(1:NTYPES,4) = NBCT2C4
         ! CY are bigger - not used
         NBCT2(5,1:NTYPES) = NBCT2CY      
         NBCT2(1:NTYPES,5) = NBCT2CY
         NBCT2(6,1:NTYPES) = NBCT2CY      
         NBCT2(1:NTYPES,6) = NBCT2CY      
         ! bases beads size 
         NBCT2(7:12,7:12) = NBCT2BASE
         ! phosphorus size       
         NBCT2(3,3) = NBCT2PP
         ! Dummy particles
         NBCT2(13,1:NTYPES) = NBCT2DUMMY
         NBCT2(1:NTYPES,13) = NBCT2DUMMY
         !effective radius adjustments
         ! P-P  and Mg++, Na+ and Cl- interactions
         NBSCORE(3,14:16) = NBRADLARGE ! P-ion interactions
         NBSCORE(14:16,3) = NBRADLARGE ! P-ion interactions
         NBSCORE(14,15:16) = NBRADLARGE ! Mg - Na,Cl interactions
         NBSCORE(15:16,14) = NBRADLARGE ! Na,Cl - Mg interactions
         NBSCORE(3,3) = NBRADLARGE   ! P-P interactions 
         NBSCORE(14,14) = NBRADLARGE ! Mg++ interactions      
      END SUBROUTINE SET_NBPARAMS
 
      !> Routine to populate all NB arrays
      SUBROUTINE SET_NBPARAMS_NEW()
         USE NAPARAMS, ONLY: SCORE_RNA
         NBCT2BASE = SCORE_RNA(25)
         NBCT2GEN = SCORE_RNA(26)
         NBCT2C4 = SCORE_RNA(27)
         NBCT2CY = SCORE_RNA(28)
         NBCT2PP = SCORE_RNA(29)

         NBRADGEN = SCORE_RNA(30)
         NBRADLARGE = SCORE_RNA(31)
         CALL SET_NBPARAMS()
      END SUBROUTINE SET_NBPARAMS_NEW
      
      !> Function to obtain 1-4 coefficient for two particles I and J
      !>
      !> @param[in] I - particle 1
      !> @param[in] J - particle 2
      !>
      !> @return The non-bonding 1-4 coefficient for the pair of particles provided.
      REAL(KIND = REAL64) FUNCTION GET_NBCOEF(I,J) RESULT(COEF)
         USE VAR_DEFS, ONLY: IAC
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: I,J !input integers
         INTEGER :: IDI, IDJ
         
         ! get the type for atoms i and j
         IDI = IAC(I)
         IDJ = IAC(J)
         COEF = NBCOEF(IDI,IDJ)      
      END FUNCTION GET_NBCOEF
   
END MODULE NBDEFS
