! This subroutine is to select whether we are to minimise in rigid body/atom coordinates
! ATOMRIGIDCOORDT is the toggling logical variable
SUBROUTINE MYLBFGS(N, M, XCOORDS, DIAGCO, EPS, MFLAG, ENERGY, ITMAX, ITDONE, RESET)
   USE COMMONS, ONLY: NATOMS, HYBRIDMINT, EPSRIGID, DEBUG, MYUNIT, GRADPROBLEMT, RMS
   USE GENRIGID, ONLY: DEGFREEDOMS, RIGIDINIT, ATOMRIGIDCOORDT, TRANSFORMCTORIGID, &
                       TRANSFORMRIGIDTOC, NRIGIDBODY
   USE PREC, ONLY: REAL64

   IMPLICIT NONE
! Arguments
   INTEGER, intent(in) :: N
   INTEGER, intent(in) :: M
   REAL(KIND = REAL64),intent(inout) :: XCOORDS(3*NATOMS)
   LOGICAL, intent(in) :: DIAGCO
   REAL(KIND = REAL64), intent(in) :: EPS
   LOGICAL, intent(inout) :: MFLAG
   REAL(KIND = REAL64), intent(out) :: ENERGY
   INTEGER, intent(in) :: ITMAX
   INTEGER, intent(inout) :: ITDONE
   LOGICAL, intent(in) :: RESET

   REAL(KIND = REAL64)      :: XRIGIDCOORDS(DEGFREEDOMS)
   REAL(KIND = REAL64)      :: EPS_TEMP


! if generalised rigid body is used, then use rigid coords, otherwise proceed as usual
! When passing rigid coords, pass zeroes from (DEGFREEDOMS+1:3*NATOMS)
! mymylbfgs is the old mylbfgs
   IF (RIGIDINIT) THEN
      ATOMRIGIDCOORDT = .FALSE.
      ! Convert to rigid body coordinates
      CALL TRANSFORMCTORIGID(XCOORDS, XRIGIDCOORDS)
      XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
      XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
      IF (HYBRIDMINT) THEN
         EPS_TEMP = EPSRIGID
      ELSE
         EPS_TEMP = EPS
      END IF
      CALL MYMYLBFGS(N, M, XCOORDS, DIAGCO, EPS_TEMP, MFLAG, ENERGY, ITMAX, ITDONE, RESET) 
      IF (DEBUG.AND.HYBRIDMINT) WRITE(MYUNIT, '(A)') ' HYBRIDMIN> Rigid body minimisation converged, switching to all-atom'
      ! Convert back to atomistic coordinates. 
      XRIGIDCOORDS(1:DEGFREEDOMS) = XCOORDS(1:DEGFREEDOMS)
      CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
      ATOMRIGIDCOORDT = .TRUE.
   END IF
   EPS_TEMP = EPS
! If we're not using rigid bodies, or we've already completed a rigid body minimisation and we're using
! hybrid minimisation, then perform a minimisation in atomistic coordinates.
   IF ((.NOT. RIGIDINIT) .OR. HYBRIDMINT) THEN
      CALL MYMYLBFGS(N, M, XCOORDS, DIAGCO, EPS_TEMP, MFLAG, ENERGY, ITMAX, ITDONE, RESET)
   END IF

   GRADPROBLEMT = .FALSE.
   IF (RMS.EQ.1.0D-100) THEN
       WRITE(MYUNIT,'(2A)') ' mylbfgs> WARNING - RMS force was set to 1.0e-100 as it was either smaller ', &
            'than this value or NaN - discarding structure'
       WRITE(MYUNIT,'(A)') '  - see debug output for further information'
       GRADPROBLEMT = .TRUE.
       MFLAG=.FALSE.
   ENDIF

END SUBROUTINE MYLBFGS


