!> @file
!> Contains RAND_ROUTINES module with facilities for random number generation

!> Module to provide random numbers
!> @brief
!> 
!> This module contains functions for normally distributed random numbers and random integers.
!> The algorithm to obtain random numbers is in random.f90.
!> The functions in random.f90, DPRAND and SPRAND are authored by N.M. Maclaren in 1992 at the University of Cambridge.
MODULE RAND_ROUTINES
   USE NUMKIND
   IMPLICIT NONE
   !> Variable to save whether we have generated already another random number, as RAND_NORMAL produces pairs of random numbers.
   LOGICAL, SAVE :: SPARET = .FALSE.
   !> Spare random number for RAND_NORMAL
   REAL(KIND = REAL64), SAVE :: SPARE = 0.0D0
   
   CONTAINS
   
   !> Subroutine to provide random normal distributed reals
   !> @brief
   !>
   !> Creation of pseudo-random normal distributed variables
   !> Polar method by Marsaglia and Bray
   !> Marsaglia, G.; Bray, T. A. (1964). 
   !> "A Convenient Method for Generating Normal Variables". SIAM Review. 6 (3): 260–264.
   !>
   !> @param[in] STDEV - standard deviation of the normal distribution to be sampled
   !> @param[in] MEAN - mean of the normal distribution to be sampled
   !> @param[out] RANDNORM - Output pseudo-random REAL64
   !>
   !> @see DPRAND
   SUBROUTINE RAND_NORMAL(STDEV, MEAN, RANDNORM)
      REAL(KIND = REAL64), INTENT(IN) :: STDEV, MEAN 
      REAL(KIND = REAL64), INTENT(OUT) :: RANDNORM
      REAL(KIND = REAL64) :: U, V, R2, S, DPRAND
      LOGICAL :: ACCEPTPAIRT
      
      ! if we have created a pair the last time, we can now use the second 
      ! number form the pair
      IF (SPARET) THEN
         SPARET = .FALSE.
         RANDNORM = SPARE * STDEV + MEAN
      ! otherwise create a new pair
      ELSE
         ACCEPTPAIRT = .FALSE.
         DO WHILE (.NOT.ACCEPTPAIRT)
            U = 2.0 * DPRAND() - 1.0
            V = 2.0 * DPRAND() - 1.0
            R2 = U * U + V * V
            IF (R2.LT.1.0D0) ACCEPTPAIRT = .TRUE.
         END DO
         S = SQRT(-2.0 * LOG(R2) / R2)
         SPARE = V * S
         SPARET = .TRUE.
         RANDNORM = U * S * STDEV + MEAN
      END IF 
   END SUBROUTINE RAND_NORMAL

   !> Reset moduel vairbale
   SUBROUTINE RESET_RANDNORM() 
      SPARET = .FALSE.
      SPARE = 0.0D0
   END SUBROUTINE RESET_RANDNORM

   !> Return random interger in range
   !> @brief
   !>
   !> Draw random integer from range LOWER to UPPER including LOWER and UPPER.
   !> We scale the random number drawn from 0...1 to the range of UPPER to LOWER+1.
   !> We then round down to the nearest integer.
   !>
   !> @param[in] LOWER - Lower end of the range
   !> @param[in] UPPER - Upper end of the range
   !> @param[out] RINT - Random integer in range
   !>
   !> @see DPRAND
   SUBROUTINE RANDINT(LOWER,UPPER,RINT)
      INTEGER, INTENT(IN)  :: LOWER, UPPER
      INTEGER, INTENT(OUT) :: RINT
      REAL(KIND = REAL64)  :: RAND, DPRAND
    
      RAND = DPRAND()
      ! scale RAND to range from LOWER to UPPER+1
      RAND = (UPPER-LOWER+1) * RAND + LOWER
      ! use internal nearest to get next integer, -1.0 goes for integer smaller than rand
      RINT = NEAREST(RAND, -1.0)
   END SUBROUTINE RANDINT

END MODULE RAND_ROUTINES
