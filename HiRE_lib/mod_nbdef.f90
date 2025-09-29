!> @file
!> Contains NBDEFS module that handles the variables for the non-bonded terms

!> The module contains the definitions and parameter values for the non-bonded terms.\n
!> It also contains the setup to populate the non-bonded scoring matrix and a function to obtain values from it.
MODULE NBDEFS
   USE PREC_HIRE
   ! Non-bonded variables
   ! NTYPES is the number of particle types, found in the .top file
   ! While the number can potentially change, we have 19 types implemented 
   ! Currently, it's taken to be :
   ! 1: O3'  2: O5'  3: P  4: R4 (RNA C4') 5: S4 (RNA C4') 6: R1 (RNA C1'), 7: S1 (DNA C1'),
   ! 8,9: G1,2  10,11: A1,2  12: U1  13: C1, 14: T1,
   ! 15: D 16: MG 17: NA  18: K 19: CL

   !> Number of different particles currently implemented
   INTEGER, PARAMETER  :: NTYPES = 19    

   !> Size paramters for large beads (P, R4, R1) 
   REAL(KIND = REAL64) :: NBCT2LARGE = 4.0D0
   !> Size paramters for medium beads (O3 and O5)
   REAL(KIND = REAL64) :: NBCT2MEDIUM = 3.6D0
   !> Size paramters for small bead (X1 and X2)
   REAL(KIND = REAL64) :: NBCT2SMALL = 3.2D0
   !> Array with bead sizes
   REAL(KIND = REAL64) :: NBCT2(NTYPES,NTYPES)

   !> Charges of particle types, these charges should never change!
   REAL(KIND = REAL64), DIMENSION(NTYPES), PARAMETER ::  &
                         CHRG = (/0.0, 0.0,-1.0, 0.0, 0.0, 0.0, &
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, &
                                  0.0, 0.0, 0.0, 2.0, 1.0, 1.0, &
                                  -1.0 /) 
                                               
   CONTAINS
      !> Routine to populate NBCT2 array
      SUBROUTINE SET_NBPARAMS()
         ! sizes for particle interactions - the larger particle is used (P and A2 for example is 4.0, while B1 and B2 is 3.2)
         ! set 4.0 as base type
         NBCT2(1:NTYPES,1:NTYPES) = NBCT2LARGE
         ! adjust sizes for base
         NBCT2(8:14,8:14) = NBCT2SMALL
         ! adjust for oxygen (O3 and O5)
         NBCT2(1:2,1:2) = NBCT2MEDIUM
         NBCT2(1:2,8:14) = NBCT2MEDIUM
         NBCT2(8:14,1:2) = NBCT2MEDIUM    
      END SUBROUTINE SET_NBPARAMS   
END MODULE NBDEFS
