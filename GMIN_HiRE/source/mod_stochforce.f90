! kr366> Module to take stochastic small forces steps
!
! We select a number of gradient components and reverse their sign and increase their magnitude.
! This chnage allows us to move away from a minimum. After some steps are taken,
! we turn of the gradient modification and complete the minimisation.
! The steps are applied during the L-BFGS minimisation.
! While the selection could be random (think SGD), we can also aim for specific motions,
! in particular by selecting the components with the samllest contributions, i.e. the softest modes.

MODULE STOCH_FORCE_STEPS
   USE PREC, ONLY: REAL64
   IMPLICIT NONE
   ! use stochastic force steps
   LOGICAL :: STOCHFORCET = .FALSE.
   ! Switch toggling whether the multiplication is applied
   LOGICAL :: STOCHFORCESTEPT = .FALSE.
   ! Array containing the multiplication factors
   REAL(KIND=REAL64), ALLOCATABLE :: FACTORS(:)
   REAL(KIND=REAL64) :: MULTIPLIER = -10.0
   ! Modulation and RMS bounds
   REAL(KIND=REAL64) :: RMSCUTOFF = 10000.0
   REAL(KIND=REAL64) :: LOWERMOD = 0.9
   REAL(KIND=REAL64) :: RAISEMOD = 1.2
   REAL(KIND=REAL64) :: JUMPMULTIPLIER = 100.0
   REAL(KIND=REAL64) :: MINMULT = 2.0
   ! Number of force modulation steps
   INTEGER :: NMAXMULTSTEPS = 0
   ! Number of components chosen
   INTEGER :: NCOMPONENTS = 0
   ! Switch to toggle between softest modes and random selections
   LOGICAL :: RANDOMFACTORST = .FALSE.
   ! Indices of chosen components (faster loop)
   INTEGER, ALLOCATABLE :: COMPSIDX(:)
   ! Switch whether we have a remaining random normal number
   LOGICAL :: SPARET = .FALSE.
   REAL(KIND=REAL64) :: SPARE
   CONTAINS

      ! subroutine to take step
      SUBROUTINE GRADMOD_STEP(COORDS)
         USE COMMONS, ONLY: NATOMS
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: COORDS(3*NATOMS)
         REAL(KIND=REAL64) :: EINIT, GRAD(3*NATOMS), RMS
         INTEGER :: J

         CALL POTENTIAL(COORDS,GRAD,EINIT,.TRUE.,.FALSE.)
         CALL GET_FACTORS(GRAD) 
         DO J=1,NMAXMULTSTEPS
            CALL POTENTIAL(COORDS,GRAD,EINIT,.TRUE.,.FALSE.)
            CALL MULTIPLY_FACTORS(GRAD)
            RMS=MAX(DSQRT(SUM(GRAD(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100)
            IF (RMS.GT.RMSCUTOFF) THEN
               MULTIPLIER = LOWERMOD*MULTIPLIER
               IF (ABS(MULTIPLIER).LT.MINMULT) THEN
                  MULTIPLIER = JUMPMULTIPLIER*MULTIPLIER
               END IF
               EXIT
            END IF
            COORDS(1:3*NATOMS) = COORDS(1:3*NATOMS) + GRAD(1:3*NATOMS)
            WRITE(*,*) " RMS: ", RMS
            IF (J.EQ.NMAXMULTSTEPS) THEN
               MULTIPLIER = RAISEMOD*MULTIPLIER
            END IF
         END DO
      END SUBROUTINE GRADMOD_STEP

      ! Subroutine to fill factor array
      SUBROUTINE GET_FACTORS(GRAD)
         USE COMMONS, ONLY: NATOMS
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: GRAD(3*NATOMS)
         ! empty any previous array and allow for a change in the number of coordinates
         ! then initialise everything with 1.0
         IF (ALLOCATED(FACTORS)) DEALLOCATE(FACTORS)
         ALLOCATE(FACTORS(3*NATOMS))
         FACTORS(1:3*NATOMS) = 1.0D0
         ! Make sure the index storage array is allocated
         IF (.NOT.(ALLOCATED(COMPSIDX))) ALLOCATE(COMPSIDX(NCOMPONENTS))
         COMPSIDX(1:NCOMPONENTS) = 0
         ! Now get the correct factors, either by random selection ...
         IF (RANDOMFACTORST) THEN
            CALL GET_RANDOM()
         ! ... or by selecting the smallest NCOMPONENTS
         ELSE
            CALL GET_SMALLEST_COMPS(GRAD)
         END IF
      END SUBROUTINE GET_FACTORS

      ! Subroutine to multiple gradient with factor
      ! We rely on the compiler to turn this simple loop into fast code,
      ! but generally this should work
      SUBROUTINE MULTIPLY_FACTORS(GRAD)
         USE COMMONS, ONLY: NATOMS
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: GRAD(3*NATOMS)
         INTEGER :: J, IDX

         DO J=1,NCOMPONENTS
            IDX = COMPSIDX(J)
            GRAD(IDX) = GRAD(IDX)*FACTORS(IDX)
         END DO
      END SUBROUTINE MULTIPLY_FACTORS

      SUBROUTINE FINISH_STOCHSTEPS()
         IF (ALLOCATED(FACTORS)) DEALLOCATE(FACTORS)
         IF (ALLOCATED(COMPSIDX)) DEALLOCATE(COMPSIDX)
      END SUBROUTINE FINISH_STOCHSTEPS

      SUBROUTINE GET_SMALLEST_COMPS(GRAD)
         USE COMMONS, ONLY: NATOMS
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: GRAD(3*NATOMS)
         INTEGER :: J

         CALL GET_LOWEST_VALS(3*NATOMS, GRAD, NCOMPONENTS, COMPSIDX)
         DO J=1,NCOMPONENTS
            FACTORS(COMPSIDX(J)) = MULTIPLIER
         ENDDO         
      END SUBROUTINE GET_SMALLEST_COMPS


      ! Select random factors
      ! First, components are randomly selected (uniform distribution)
      ! Then, for each a random factor is assigned drawn from a normal distribution
      SUBROUTINE GET_RANDOM()
         IMPLICIT NONE
         INTEGER :: J
         REAL(KIND=REAL64) :: RANDOUT
         REAL(KIND=REAL64) :: STDEV
         REAL(KIND=REAL64) :: MEAN
         ! get n random indices
         CALL SELECT_NRANDS()
         ! iterate and get a random facto
         STDEV = ABS(0.15*MULTIPLIER)
         MEAN = ABS(MULTIPLIER)
         DO J=1,NCOMPONENTS
            ! get random number from normal distribution
            CALL RAND_NORMAL(STDEV,MEAN,RANDOUT)
            FACTORS(COMPSIDX(J)) = -RANDOUT
         ENDDO
      END SUBROUTINE GET_RANDOM

      ! select NCOMPONENTS random indices (we throw away duplicates, so NOPT >> NCOMPONENTS is desirable)
      SUBROUTINE SELECT_NRANDS()
         USE COMMONS, ONLY: NATOMS
         IMPLICIT NONE
         LOGICAL :: PICKED(3*NATOMS)
         INTEGER :: NPICKED, RANDIDX

         ! No value has been chosen at this point
         PICKED(1:3*NATOMS) = .FALSE.
         NPICKED = 0
         DO WHILE (NPICKED.LT.NCOMPONENTS)
            ! get a random integer between 1 and nopt
            CALL RANDINT(3*NATOMS,RANDIDX)
            ! if this number was already picked, we simply cycle through
            IF (PICKED(RANDIDX)) THEN
               CYCLE
            ! otherwise, we add the new entry to our list
            ELSE
               NPICKED = NPICKED + 1
               COMPSIDX(NPICKED) = RANDIDX
               PICKED(RANDIDX) = .TRUE.
            ENDIF
         END DO
      END SUBROUTINE SELECT_NRANDS

      ! Draw random integer between 1 and UPPER
      SUBROUTINE RANDINT(UPPER,RINT)
         INTEGER, INTENT(IN)  :: UPPER
         INTEGER, INTENT(OUT) :: RINT
         REAL(KIND = REAL64)  :: RAND, DPRAND
         
         RAND = DPRAND()
         ! scale RAND to range from LOWER to UPPER+1
         RAND = UPPER * RAND + 1
         ! use internal nearest to get next integer, -1.0 goes for integer smaller than rand
         RINT = NEAREST(RAND, -1.0)
      END SUBROUTINE RANDINT

      ! Creation of pseudo-random normal distributed variables
      ! Polar method by Marsaglia and Bray
      ! Marsaglia, G.; Bray, T. A. (1964). 
      ! "A Convenient Method for Generating Normal Variables". SIAM Review. 6 (3): 260–264.
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

      SUBROUTINE GET_LOWEST_VALS(NDIM, ARR, K, OUTARR)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NDIM
         INTEGER, INTENT(IN) :: K
         REAL(KIND=REAL64), INTENT(IN) :: ARR(NDIM) 
         INTEGER, INTENT(OUT) :: OUTARR(K)
         REAL(KIND=REAL64) :: KVAL
         INTEGER :: J, NDUMMY
         REAL(KIND=REAL64), PARAMETER :: EPS = 1.0D-11
         KVAL = SELECT_KTH(NDIM, ARR, K)
         NDUMMY = 0
         OUTARR(1:K) = 0
         DO J=1,NDIM
            IF (ARR(J).LE.KVAL) THEN
               NDUMMY = NDUMMY + 1
               OUTARR(NDUMMY) = J
            END IF
            IF (NDUMMY.EQ.K) EXIT
         END DO

         IF (NDUMMY.LT.K) THEN
            DO J=1,NDIM
                IF (ABS(ARR(J)-KVAL).LT.EPS) THEN
                   WRITE(*,*) "J, NDUMMY:", J, NDUMMY
                   WRITE(*,*) OUTARR
                   NDUMMY = NDUMMY + 1
                   OUTARR(NDUMMY) = J
                END IF
                IF (NDUMMY.EQ.K) EXIT
             END DO
         ENDIF 
      END SUBROUTINE GET_LOWEST_VALS

      REAL(KIND=REAL64) FUNCTION SELECT_KTH(NDIM,INARR,K)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NDIM                  ! array dimension
         REAL(KIND=REAL64), INTENT(IN) :: INARR(NDIM) ! array 
         INTEGER, INTENT(IN) :: K                     ! Kth element
         INTEGER :: L1, R1                            ! outer fingers
         INTEGER :: L2, R2                            ! inner fingers
         REAL(KIND=REAL64) :: EL1, EL2                ! elements tested currently
         REAL(KIND=REAL64) :: ARR(NDIM)               ! working array
         REAL(KIND=REAL64), PARAMETER :: EPS = 1.0D-12
         INTEGER :: NDUMMY
         NDUMMY = 0
         ARR(1:NDIM) = INARR(1:NDIM)
         ! set up bounds
         L1 = 1
         R1 = NDIM
         ! WRITE(*,*) "L1 and R1: ", L1, R1
         ! run until the element is clamped from both sites 
         DO WHILE (L1 .LT. R1)
            EL1 = ARR(K)   ! for a sorted array, this would be the ideal case
            ! set inner teeth for scan
            L2 = L1
            R2 = R1
            ! WRITE(*,*) "L2 and R2: ", L2, R2
            DO WHILE (L2.LE.R2)
               ! parse elements less than EL1 and raise left inner tooth
               DO WHILE (ARR(L2).LT.EL1)
                  L2 = L2 + 1
               END DO
               ! parse elements larger than EL1 and lower right inner tooth
               DO WHILE (EL1.LT.ARR(R2))
                  R2 = R2 -1
               END DO
               ! check where the teeth are at
               IF (L2.LT.R2) THEN
                  ! neither tooth has reached the pivot, switch the elements
                  EL2 = ARR(L2)
                  ARR(L2) = ARR(R2)
                  ARR(R2) = EL2
               ELSE IF (L2.EQ.R2) THEN
                  ! push the teeth one further to leave do while loop
                  L2 = L2 + 1
                  R2 = R2 - 1
               END IF
               ! WRITE(*,*) "EL1 and EL2: ", EL1, EL2
               IF (EL1.EQ.EL2) THEN
                  ARR(R2) = ARR(R2) + ((-1)**NDUMMY)*EPS
                  NDUMMY = NDUMMY + 1
               END IF
            END DO
            IF (R2.LT.K) THEN
               L1 = L2
            END IF
            IF (K.LT.L2) THEN
               R1 = R2
            END IF
            ! WRITE(*,*) "L1 and R1 at end of outer loop: ", L1, R1
         END DO
         SELECT_KTH = ARR(K)
         ! WRITE(*,*) ARR
      END FUNCTION SELECT_KTH

END MODULE STOCH_FORCE_STEPS
