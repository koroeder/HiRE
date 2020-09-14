MODULE MOVES
! This module is for the magical new way of doing custom movesets and also fixing our
! horrible old way of doing moves with 15 different implementations of Cartesian steps.
! Please code nicely in here...GOTOs will be removed, as will unclear variable names.
! Consistent indentation is mandatory!
USE COMMONS
USE PREC
IMPLICIT NONE

CONTAINS

SUBROUTINE CARTESIAN_SPHERE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random spherically symmetric displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved

! The VEC3 module (vec3.f90) contains helper functions for handling vectors and matricies
   USE VEC3
! The SANITY module contains sanity check functions
   USE SANITY
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   REAL(KIND = REAL64)                              :: DPRAND
   REAL(KIND = REAL64), INTENT(IN)                  :: MAX_STEP
   REAL(KIND = REAL64), DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SPHERE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
      XYZ(3*I-2:3*I)=XYZ(3*I-2:3*I)+VEC_RANDOM()*(DPRAND()**(1.0D0/3.0D0))*MAX_STEP
   ENDDO

END SUBROUTINE CARTESIAN_SPHERE

SUBROUTINE CARTESIAN_SIMPLE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved

! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   REAL(KIND = REAL64)                              :: DPRAND
   REAL(KIND = REAL64)                              :: RANDOMX
   REAL(KIND = REAL64)                              :: RANDOMY
   REAL(KIND = REAL64)                              :: RANDOMZ
   REAL(KIND = REAL64), INTENT(IN)                  :: MAX_STEP
   REAL(KIND = REAL64), DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SIMPLE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
! Draw a random number between -1 and 1 for each coordinate
      RANDOMX=(DPRAND()-0.5D0)*2.0D0
      RANDOMY=(DPRAND()-0.5D0)*2.0D0
      RANDOMZ=(DPRAND()-0.5D0)*2.0D0
! Displace each coordinate
      XYZ(3*I-2)=XYZ(3*I-2)+MAX_STEP*RANDOMX
      XYZ(3*I-1)=XYZ(3*I-1)+MAX_STEP*RANDOMY
      XYZ(3*I  )=XYZ(3*I  )+MAX_STEP*RANDOMZ
   ENDDO

END SUBROUTINE CARTESIAN_SIMPLE

SUBROUTINE ROTATION_ABOUT_AXIS(XYZ, VECTOR_START_XYZ, &
                               VECTOR_END_XYZ, ANGLE_DEGS, ATOM_LIST)
!
! Rotate the coordinates of the atoms in ATOM_LIST about the line from VECTOR_START_XYZ
! to VECTOR_END_XYZ through ANGLE degrees.
!
! Arguments
! ---------
!
! Required:
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! VECTOR_START_XYZ(in): start of the line about which to rotate, in Cartesian coords
! VECTOR_END_XYZ(in): end of the line about which to rotate, in Cartesian coords
! ANGLE_DEGS(in): angle through which to rotate, in degrees
!
! Optional:
! ATOM_LIST(in): list of atoms to be rotated
!
! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
! Arguments
   REAL(KIND = REAL64), DIMENSION(:), INTENT(INOUT)   :: XYZ
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)     :: ATOM_LIST
   REAL(KIND = REAL64), DIMENSION(3), INTENT(IN)      :: VECTOR_START_XYZ
   REAL(KIND = REAL64), DIMENSION(3), INTENT(IN)      :: VECTOR_END_XYZ
   REAL(KIND = REAL64), INTENT(IN)                    :: ANGLE_DEGS
! Constants
! Variables
   REAL(KIND = REAL64)                                :: PI
   REAL(KIND = REAL64)                                :: DEGS_OVER_RADS
   REAL(KIND = REAL64)                                :: COS_THETA, SIN_THETA
   REAL(KIND = REAL64)                                :: A, B, C
   REAL(KIND = REAL64)                                :: D, E, F
   REAL(KIND = REAL64)                                :: U, V, W
   REAL(KIND = REAL64)                                :: X, Y, Z
   REAL(KIND = REAL64)                                :: VECTOR_MAG
   INTEGER                                         :: NUM_ATOMS
   INTEGER                                         :: I
   LOGICAL, ALLOCATABLE, DIMENSION(:)              :: ATOM_MASK
   LOGICAL                                         :: TEST
! Function declarations
   REAL(KIND = REAL64)                                :: DNRM2

! Define PI using ATAN and the conversion factor between degrees and radians.
! x degrees = x * DEGS_OVER_RADS radians
   PI = 4.0D0 * ATAN(1.0D0)
   DEGS_OVER_RADS = PI / 180.0D0

! Calculate cos and sin of the angle.
   SIN_THETA = SIN(ANGLE_DEGS * DEGS_OVER_RADS)
   COS_THETA = COS(ANGLE_DEGS * DEGS_OVER_RADS)

! Assign variables according to the functional form described on:
! http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html

   A = VECTOR_START_XYZ(1)
   B = VECTOR_START_XYZ(2)
   C = VECTOR_START_XYZ(3)
   
   D = VECTOR_END_XYZ(1)
   E = VECTOR_END_XYZ(2)
   F = VECTOR_END_XYZ(3)

! DNRM2(INTEGER N, REAL X(N), INTEGER INCX)
   VECTOR_MAG = DNRM2(3, VECTOR_END_XYZ - VECTOR_START_XYZ, 1)

   U = (D - A) / VECTOR_MAG
   V = (E - B) / VECTOR_MAG
   W = (F - C) / VECTOR_MAG

! Work out how many atoms the coordinates passed describe
   NUM_ATOMS = SIZE(XYZ) / 3 
! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to ROTATION_ABOUT_AXIS'
   ENDIF

! Convert ATOM_LIST (a list of atoms to be rotated) into ATOM_MASK, which can
! applied in a WHERE loop.
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Loop through and apply the formula if the atom described is in ATOM_LIST
! N.B. I've rearranged the formula so that the variables cycle, since I think it
! makes checking a bit easier.
   DO I = 1, NUM_ATOMS
      IF (.NOT. ATOM_MASK(I)) CYCLE
      X = XYZ(3 * I - 2)
      Y = XYZ(3 * I - 1)
      Z = XYZ(3 * I    )
      XYZ(3 * I - 2) = (A * (V**2 + W**2) - U * (B*V + C*W - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       X * COS_THETA + &
                       (-C*V + B*W + V*Z - W*Y) * SIN_THETA
      XYZ(3 * I - 1) = (B * (W**2 + U**2) - V * (C*W + A*U - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Y * COS_THETA + &
                       (-A*W + C*U + W*X - U*Z) * SIN_THETA
      XYZ(3 * I    ) = (C * (U**2 + V**2) - W * (A*U + B*V - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Z * COS_THETA + &
                       (-B*U + A*V + U*Y - V*X) * SIN_THETA
   END DO

END SUBROUTINE ROTATION_ABOUT_AXIS

   SUBROUTINE PULL_MOVE(X)
      REAL(KIND=REAL64) , INTENT(INOUT) :: X(3*NATOMS)
      REAL(KIND=REAL64) :: EDUMMY, RANDOM, DPRAND
      INTEGER :: U, V, END1, END2, ITER=0
      LOGICAL :: CFLAG    
      !pull molecule and then release
      !pick nucleotides in the first and last third of molecule
      END1 = NNUCL/3
      END2 = 2*NNUCL/3
      RANDOM = DPRAND()
      U = FLOOR(END1*RANDOM+1)
      RANDOM = DPRAND()
      V = FLOOR((NNUCL-END2+1)*RANDOM+1)
      PULLT = .TRUE.
      PATOM1 = LIST_NUCL(U)%FATOM
      PATOM2 = LIST_NUCL(V)%FATOM
      PFORCE = PULLMF
      WRITE(MYUNIT, '(A,I8,A,I8)') " moves> Apply pulling force for nucleotides ", U , " and " , V
      WRITE(MYUNIT, '(A,F12.3)') " moves> Pulling force is: ", PFORCE
      CALL MYLBFGS(3*NATOMS,MUPDATE,X,.FALSE.,BQMAX,CFLAG,EDUMMY,MAXIT,ITER,.TRUE.)     
      PULLT = .FALSE. 
      PATOM1 = 0
      PATOM2 = 0
      RETURN
   END SUBROUTINE PULL_MOVE

   SUBROUTINE TWIST_MOVE(X)
      USE TWIST_MOD
      REAL(KIND=REAL64) , INTENT(INOUT) :: X(3*NATOMS)
      REAL(KIND=REAL64) :: EDUMMY, RANDOM, DPRAND
      INTEGER :: U, V, END1, END2, ITER=0
      LOGICAL :: CFLAG          
      !Twist molecule and then release
      END1 = NNUCL/3
      END2 = 2*NNUCL/3     
      RANDOM = DPRAND()-0.5 !get random number between -0.5 and 0.5
      U = END1 + FLOOR(5*RANDOM+0.5)
      RANDOM = DPRAND()-0.5 !get random number between -0.5 and 0.5
      V = END2 + FLOOR(5*RANDOM+0.5) 
      !use the twist module for this
      NTWISTGROUPS=1
      ALLOCATE(TWIST_K(1))
      ALLOCATE(TWIST_THETA0(1))
      ALLOCATE(TWIST_ATOMS(4,1))
      TWIST_ATOMS(1,1) = LIST_NUCL(1)%FATOM
      TWIST_ATOMS(2,1) = LIST_NUCL(U)%FATOM
      TWIST_ATOMS(3,1) = LIST_NUCL(V)%FATOM
      TWIST_ATOMS(4,1) = LIST_NUCL(NNUCL)%FATOM
      TWIST_K(1) = TWISTMF
      RANDOM = DPRAND()
      TWIST_THETA0(1) = 2*PI*(RANDOM-0.5) !random equilibrium angle
      TWISTT = .TRUE.
      WRITE(MYUNIT, '(A,I8,A,I8,A,I8)') " moves> Apply twist for nucleotides 1, ", U , ", " , V, " and ", NNUCL     
      WRITE(MYUNIT, '(A,F12.3,A,F8.4)') " moves> Applied force: ", TWIST_K(1), ", eq. angle: ", TWIST_THETA0(1)
      CALL MYLBFGS(3*NATOMS,MUPDATE,X,.FALSE.,BQMAX,CFLAG,EDUMMY,MAXIT,ITER,.TRUE.)    
      TWISTT = .FALSE.
      DEALLOCATE(TWIST_K)
      DEALLOCATE(TWIST_THETA0)
      DEALLOCATE(TWIST_ATOMS)
      RETURN
   END SUBROUTINE TWIST_MOVE

   SUBROUTINE HARMONIC_MOVE(X,MOVEDT)
      REAL(KIND=REAL64) , INTENT(INOUT) :: X(3*NATOMS)
      LOGICAL, INTENT(OUT) :: MOVEDT
      REAL(KIND=REAL64) :: CDUMMY(3*NATOMS), EDUMMY, RANDOM, DPRAND
      LOGICAL :: INBP(NNUCL), CFLAG
      INTEGER :: I1, I2, I3, STARTSTOP(NNUCL,2), STARTID, ENDID, NSTRETCH, U, V, N1, N2, ITER=0
      !add spring between two free nucleotides and quench, then remove spring and pass back
      MOVEDT=.FALSE.
      ! 1. find current base pairs
      INBP(:) = .FALSE.
      DO I1=1,NNUCL-1
         DO I2=2,NNUCL
            IF (BP_CURR(I1,I2)) THEN
               INBP(I1) = .TRUE.
               INBP(I2) = .TRUE.
            ENDIF
         ENDDO
      ENDDO
      ! 2. now determine stretches without base pairs (a stretch must be 4 nucleotides or longer)
      STARTID=-1
      ENDID=-1
      NSTRETCH=0

      DO I3=1,NNUCL
         IF (.NOT.INBP(I3)) THEN
            IF (STARTID.EQ.-1) THEN !if start isn't set, this is the first free nucleotide
               STARTID=I3
            ELSE                    !otherwise update ENDID
               ENDID=I3
            ENDIF
         ELSE
            IF ((STARTID.NE.-1).AND.(ENDID.NE.-1)) THEN ! only if both are set we might have a stretch
               IF (DEBUG) WRITE(MYUNIT,*) " moves> free stretch between ", STARTID, " and ", ENDID
               IF ((ENDID-STARTID).GT.3) THEN           ! apply minimum length condition
                  NSTRETCH = NSTRETCH + 1 
                  STARTSTOP(NSTRETCH,1) = STARTID
                  STARTSTOP(NSTRETCH,2) = ENDID
               ENDIF
            ENDIF
            STARTID = -1
            ENDID = -1
         ENDIF
      ENDDO
      !The last stretch wouldn't be closed, so we need to
      IF (ENDID.NE.-1) THEN
         IF ((ENDID-STARTID).GT.3) THEN
            NSTRETCH = NSTRETCH + 1
            STARTSTOP(NSTRETCH,1) = STARTID
            STARTSTOP(NSTRETCH,2) = ENDID
         ENDIF
      ENDIF

      IF (DEBUG) WRITE(MYUNIT, '(A,I8)') " moves> Number of unpaired stretches found: ", NSTRETCH
      ! test whether we have any stretches at all
      IF (NSTRETCH.EQ.0) RETURN
      ! select two nuleotides U and V to be moved with spring
      ! if we have one stretch only try to connect two nucleotides on either side of the loop
      IF (NSTRETCH.EQ.1) THEN
         STARTID = STARTSTOP(1,1)
         ENDID = FLOOR(STARTSTOP(1,2)/2.0 + STARTSTOP(1,1)/2.0) - 1
         IF (STARTID.GE.ENDID) THEN
            U = STARTID
         ELSE
            RANDOM = DPRAND()
            U = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
         ENDIF
         STARTID = ENDID + 2
         ENDID = STARTSTOP(1,2)
         IF (STARTID.GE.ENDID) THEN
            V = ENDID
         ELSE
            RANDOM = DPRAND()
            V = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
         ENDIF
      ! if we have two select a random nucleotide in both and connect
      ELSE IF (NSTRETCH.EQ.2) THEN
         ! select a random number from 0 to 1:
         RANDOM = DPRAND()
         ! now generate a random integer U between STARTID and ENDID for nucleotide 1
         STARTID = STARTSTOP(1,1)
         ENDID = STARTSTOP(1,2)
         U = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
         ! repeat for nucleotide two using V
         STARTID = STARTSTOP(2,1)
         ENDID = STARTSTOP(2,2)
         V = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
      ! if we have more than two regions, select two at random (but make sure we select two!)
      ELSE
         ! first select random stretches
         RANDOM=DPRAND()
         N1 = FLOOR(NSTRETCH*RANDOM+1)
         N2 = N1
         DO WHILE (N1.EQ.N2)
             RANDOM=DPRAND()
             N2 = FLOOR(NSTRETCH*RANDOM+1)
         END DO
         ! now find U and V
         STARTID = STARTSTOP(N1,1)
         ENDID = STARTSTOP(N1,2)
         U = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
         ! repeat for nucleotide two using V
         STARTID = STARTSTOP(N2,1)
         ENDID = STARTSTOP(N2,2)
         V = FLOOR((ENDID-STARTID+1)*RANDOM+STARTID)
      ENDIF
      HARMATOMS(1) = LIST_NUCL(U)%LATOM
      HARMATOMS(2) = LIST_NUCL(V)%LATOM               
      
      ! we add a harmonic potential between the atoms, and minimise
      WRITE(MYUNIT, '(A,I8,A,I8)') " moves> Add spring between nucleotides ", U , " and " , V
      HARMONICPOT = .TRUE.
      CDUMMY(:) = X(:)
      CALL MYLBFGS(3*NATOMS,MUPDATE,CDUMMY,.FALSE.,BQMAX,CFLAG,EDUMMY,MAXIT,ITER,.TRUE.)
      IF (CFLAG) THEN !if converged do some tests
         WRITE(MYUNIT, '(A)') " moves> Converged minimisation, remove spring"
         X(:) = CDUMMY(:)
         MOVEDT=.TRUE.
      ELSE
         ! attept minimisation again with more steps
         CALL MYLBFGS(3*NATOMS,MUPDATE,CDUMMY,.FALSE.,BQMAX,CFLAG,EDUMMY,2*MAXIT,ITER,.FALSE.)
         IF (CFLAG) THEN
            WRITE(MYUNIT, '(A)') " moves> Converged minimisation, remove spring"
            X(:) = CDUMMY(:)
            MOVEDT=.TRUE.
         ELSE
            WRITE(MYUNIT, '(A)') " moves> Could not converge structure, no move attempted"
         ENDIF
      ENDIF
      HARMONICPOT = .FALSE.
      RETURN
   END SUBROUTINE HARMONIC_MOVE

   SUBROUTINE ADAPT_FORCE(DOTWIST,DOPULL,ATEST)
      LOGICAL, INTENT(IN) :: DOTWIST, DOPULL, ATEST

      ! If move was accepted, increase force for next time to force larger change
      ! If move was rejected, decrease force
      ! Stay within the upper and lower limits 
      IF (DOTWIST) THEN
         IF (ATEST) THEN
            TWISTMF = TWISTMF*TADAPTSCALE
            IF (TWISTMF.GT.TUPPERF) TWISTMF=TUPPERF
            WRITE(MYUNIT, '(A,F12.3)') " adaptf> Increased twisting force to ", TWISTMF
         ELSE
            TWISTMF = TWISTMF/TADAPTSCALE
            IF (TWISTMF.LT.TLOWERF) TWISTMF=TLOWERF
            WRITE(MYUNIT, '(A,F12.3)') " adaptf> Decreased twisting force to ", TWISTMF
         ENDIF
      ENDIF
      IF (DOPULL) THEN
         IF (ATEST) THEN
            PULLMF = PULLMF*PADAPTSCALE
            IF (PULLMF.GT.PUPPERF) PULLMF=PUPPERF
            WRITE(MYUNIT, '(A,F12.3)') " adaptf> Increased pulling force to ", PULLMF
         ELSE
            PULLMF = PULLMF/PADAPTSCALE
            IF (PULLMF.LT.PLOWERF) PULLMF=PLOWERF
            WRITE(MYUNIT, '(A,F12.3)') " adaptf> Decreased pulling force to ", PULLMF
         ENDIF
      ENDIF           
      RETURN
   END SUBROUTINE ADAPT_FORCE

END MODULE MOVES
