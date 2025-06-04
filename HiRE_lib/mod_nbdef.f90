!> @file
!> Contains NBDEFS module that handles the variables for the non-bonded terms

!> The module contains the definitions and parameter values for the non-bonded terms.\n
!> It also contains the setup to populate the non-bonded scoring matrix and a function to obtain values from it.
MODULE NBDEFS
   USE PREC_HIRE
   ! Non-bonded variables
   ! NTYPES is the number of particle types, found in the .top file
   ! While the number can potentially change, we have 17 types implemented 
   ! Currently, it's taken to be :
   ! 1: O3'  2: O5'  3: P  4: R4 (RNA C4') 5: S4 (RNA C4') 6: R1 (RNA C1'), 7: S1 (DNA C1'),
   ! 8,9: G1,2  10,11: A1,2  12: U1  13: C1, 14: T1,
   ! 15: D 16: MG 17: NA  18: CL
   !> Number of different particles currently implemented
   INTEGER, PARAMETER  :: NTYPES = 18              
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
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                  0.0, 0.0, 0.0, 2.0, 1.0,-1.0 /) 
                                               
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
         NBCT2(5,1:NTYPES) = NBCT2C4
         NBCT2(1:NTYPES,5) = NBCT2C4       
         ! CY are bigger - not used
         NBCT2(6,1:NTYPES) = NBCT2CY      
         NBCT2(1:NTYPES,6) = NBCT2CY
         NBCT2(7,1:NTYPES) = NBCT2CY      
         NBCT2(1:NTYPES,7) = NBCT2CY      
         ! bases beads size 
         NBCT2(8:14,8:14) = NBCT2BASE
         ! phosphorus size       
         NBCT2(3,3) = NBCT2PP
         ! Dummy particles
         NBCT2(15,1:NTYPES) = NBCT2DUMMY
         NBCT2(1:NTYPES,15) = NBCT2DUMMY
         !effective radius adjustments
         ! P-P  and Mg++, Na+ and Cl- interactions
         NBSCORE(3,16:18) = NBRADLARGE ! P-ion interactions
         NBSCORE(16:18,3) = NBRADLARGE ! P-ion interactions
         NBSCORE(16,17:18) = NBRADLARGE ! Mg - Na,Cl interactions
         NBSCORE(17:18,16) = NBRADLARGE ! Na,Cl - Mg interactions
         NBSCORE(3,3) = NBRADLARGE   ! P-P interactions 
         NBSCORE(16,16) = NBRADLARGE ! Mg++ interactions      
      END SUBROUTINE SET_NBPARAMS
 
      !> Routine to populate all NB arrays
      SUBROUTINE SET_NBPARAMS_NEW()
         USE NAPARAMS, ONLY: SCORE_RNA
         NBCT2BASE = SCORE_RNA(44)
         NBCT2GEN = SCORE_RNA(45)
         NBCT2C4 = SCORE_RNA(46)
         NBCT2CY = SCORE_RNA(47)
         NBCT2PP = SCORE_RNA(48)

         NBRADGEN = SCORE_RNA(49)
         NBRADLARGE = SCORE_RNA(50)
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
