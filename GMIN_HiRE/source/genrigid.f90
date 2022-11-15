MODULE GENRIGID

      USE PREC
      IMPLICIT NONE
      INTEGER :: NRIGIDBODY, DEGFREEDOMS, MAXSITE, NRELAXRIGIDR, NRELAXRIGIDA
      INTEGER :: XNATOMS
      INTEGER, ALLOCATABLE :: NSITEPERBODY(:), REFVECTOR(:), RIGIDSINGLES(:)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: RIGIDGROUPS
      REAL(KIND = REAL64), ALLOCATABLE :: RIGIDCOORDS(:)
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:,:,:) :: SITESRIGIDBODY
      REAL(KIND = REAL64), ALLOCATABLE :: GR_WEIGHTS(:) ! weights for com calculation, e.g. masses
      LOGICAL :: RIGIDINIT, ATOMRIGIDCOORDT, RELAXRIGIDT
      LOGICAL :: GENRIGIDT

      REAL(KIND = REAL64), ALLOCATABLE :: IINVERSE(:,:,:)

      LOGICAL :: RIGIDOPTIMROTAT, FREEZERIGIDBODYT
      REAL(KIND = REAL64) :: OPTIMROTAVALUES(3)
      LOGICAL :: AACONVERGENCET
      INTEGER, ALLOCATABLE :: LRBNPERMGROUP(:), LRBNPERMSIZE(:,:), LRBPERMGROUP(:,:), LRBNSETS(:,:), LRBSETS(:,:,:)


!   jdf43:  RIGIDISRIGID = logical array, size NATOMS, TRUE if atom is part of
!           RB
      LOGICAL, ALLOCATABLE :: RIGIDISRIGID(:)
      INTEGER, ALLOCATABLE :: RB_BY_ATOM(:) ! sn402: records which RB an atom belongs to
!   jdf43:  FROZENRIGIDBODY = logical array, size NATOMS, TRUE if RB is frozen
      LOGICAL, ALLOCATABLE :: FROZENRIGIDBODY(:)

!-----------------------------------------------------------------------------------!
! NRIGIDBODY  = number of rigid bodies
! DEGFREEDOMS = number of degrees of freedom = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS
! MAXSITE     = maximum number of sites in a rigid body
! NRELAXRIGIDR = rigid body minimisation for this number of steps
! NRELAXRIGIDA = atom minimisation for this number of steps
! NSITEPERBODY= number of rigid body sites, no need to be the same for all bodies
! REFVECTOR   = reference vector for the atomistic to rigic coordinate transformation
! RIGIDSINGLES= list of atoms not in rigid bodies
! RIGIDGROUPS = list of atoms in rigid bodies, need a file called rbodyconfig
! RIGIDCOORDS = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS coordinates
! SITESRIGIDBODY = coordinates of the rigid body sites
! RIGIDINIT   = logical variable for generalised rigid body
! ATOMRIGIDCOORDT, .TRUE. = atom coords active, .FALSE. = rigid coords active, used in mylbfgs & potential
! GENRIGIDT = generalised rigid body takestep taken if .TRUE.
!-----------------------------------------------------------------------------------!


CONTAINS

!-------------------------------------------
! Initializes basic structures
! has to be the first call to GENRIGID function in order to setup basic structures.
! After that, the array which defines the sites can be filled. Then GENRIGID_INITIALIZE
! completes the initialization of rigid bodies.
!-------------------------------------------
SUBROUTINE GENRIGID_ALLOCATE(NEW_NRIGIDBODY,NEW_MAXSITE)
  USE COMMONS, only: NATOMS, ATMASS
  USE HIRE_INTERFACE, ONLY : PASS_HIRE_MASSES

  IMPLICIT NONE
  INTEGER, intent(in) :: NEW_NRIGIDBODY, NEW_MAXSITE
  
  NRIGIDBODY = NEW_NRIGIDBODY
  MAXSITE = NEW_MAXSITE

  ! Allocate NSITEPERBODY
  IF (ALLOCATED(NSITEPERBODY)) DEALLOCATE(NSITEPERBODY)
  IF (ALLOCATED(SITESRIGIDBODY)) DEALLOCATE(SITESRIGIDBODY)
  IF (ALLOCATED(RIGIDGROUPS)) DEALLOCATE(RIGIDGROUPS)
  IF (ALLOCATED(REFVECTOR)) DEALLOCATE(REFVECTOR)
  IF (ALLOCATED(GR_WEIGHTS)) DEALLOCATE(GR_WEIGHTS)
  IF (ALLOCATED(IINVERSE)) DEALLOCATE(IINVERSE)
  ALLOCATE (NSITEPERBODY(NRIGIDBODY))
  ALLOCATE (SITESRIGIDBODY(MAXSITE,3,NRIGIDBODY))
  ALLOCATE (RIGIDGROUPS(MAXSITE,NRIGIDBODY))
  ALLOCATE (REFVECTOR(NRIGIDBODY))
  ALLOCATE (GR_WEIGHTS(NATOMS))
  ALLOCATE (IINVERSE(NRIGIDBODY,3,3))
  IF (ALLOCATED(RIGIDISRIGID)) DEALLOCATE(RIGIDISRIGID)
  IF (ALLOCATED(RB_BY_ATOM)) DEALLOCATE(RB_BY_ATOM)
  ALLOCATE (RIGIDISRIGID(NATOMS))
  ALLOCATE (RB_BY_ATOM(NATOMS))
  RIGIDISRIGID=.FALSE.
  RB_BY_ATOM(:) = -1

  CALL PASS_HIRE_MASSES(NATOMS,ATMASS)
  GR_WEIGHTS(1:NATOMS)=ATMASS(1:NATOMS)

END SUBROUTINE

!-------------------------------------------
! Setup rigid body stuff after site definitions are done
!-------------------------------------------
SUBROUTINE GENRIGID_INITIALISE(INICOORDS)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2, J3, DUMMY
  REAL(KIND = REAL64) :: XMASS, YMASS, ZMASS, PNORM, MASS
  LOGICAL :: SATOMT, RTEST
  REAL(KIND = REAL64) INICOORDS(3*NATOMS)
  
  INTEGER          :: INFO
  INTEGER, PARAMETER :: LWORK = 1000000 ! the dimension is set arbitrarily
  INTEGER :: I, J
  REAL(KIND = REAL64) :: DR(3), KBLOCK(3,3), KBEGNV(3)
  REAL(KIND = REAL64) :: WORK(LWORK)

  REAL(KIND = REAL64) :: DET
  ! SITESRIGIDBODY has to be initialised to zeroes, for some reason I have a couple of NaNs in the array (gfortran)
  SITESRIGIDBODY(:,:,:)=0.0D0

  ! initialize coordinates for rigid bodies
  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        SITESRIGIDBODY(J2,:,J1) = INICOORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO

  ! determine number of degrees of freedom
  DEGFREEDOMS = 0
  DO J1 = 1, NRIGIDBODY
     DEGFREEDOMS = DEGFREEDOMS + NSITEPERBODY(J1)
  ENDDO
  DEGFREEDOMS = 6 * NRIGIDBODY + 3 * (NATOMS - DEGFREEDOMS)

  ! Allocate further data 
  IF (ALLOCATED(RIGIDSINGLES)) DEALLOCATE(RIGIDSINGLES)
  IF (ALLOCATED(RIGIDCOORDS)) DEALLOCATE(RIGIDCOORDS)
  ALLOCATE (RIGIDSINGLES((DEGFREEDOMS/3 - 2 * NRIGIDBODY)))
  ALLOCATE (RIGIDCOORDS(DEGFREEDOMS))

  DUMMY = 0
  DO J1 = 1, NATOMS
     SATOMT = .TRUE.
     DO J2 = 1, NRIGIDBODY
        DO J3 = 1, NSITEPERBODY(J2)
           IF (J1 == RIGIDGROUPS(J3,J2)) SATOMT = .FALSE.
        ENDDO
     ENDDO
     IF (SATOMT) THEN
        IF (DUMMY.EQ.DEGFREEDOMS/3-2*NRIGIDBODY) THEN
           WRITE(*,*) "genrigid> Error. More free atoms than expected."
           WRITE(*,*) "Likely problem with rbodyconfig."
           STOP
        ENDIF
        DUMMY = DUMMY + 1
        RIGIDSINGLES(DUMMY) = J1
     ENDIF
  ENDDO
 
  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0d0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + SITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + SITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + SITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS  = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        SITESRIGIDBODY(J2,1,J1) = SITESRIGIDBODY(J2,1,J1) - XMASS
        SITESRIGIDBODY(J2,2,J1) = SITESRIGIDBODY(J2,2,J1) - YMASS
        SITESRIGIDBODY(J2,3,J1) = SITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO

  ENDDO

  DO J1 = 1, NRIGIDBODY
     IINVERSE(J1,:,:) = 0.0D0
  END DO

  IF (AACONVERGENCET .EQV. .TRUE.) THEN
     DO J1 = 1, NRIGIDBODY
        KBLOCK(:,:) = 0.0D0
        DO J2 = 1, NSITEPERBODY(J1)        
           DR(:)  = SITESRIGIDBODY(J2,:,J1)
           DO I = 1, 3
              KBLOCK(I,I) = KBLOCK(I,I) + (DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
        CALL DSYEV('V','U',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
        CALL RBDET(KBLOCK, DET)
        IF (DET < 0.0D0) THEN
           KBLOCK(:,3) = -KBLOCK(:,3)
           CALL RBDET(KBLOCK, DET)
           IF (DET < 0.0D0) THEN
              PRINT *, "GENRIGID> BAD ALIGNMENT", J1
              STOP
           ENDIF
        ENDIF
        KBLOCK = TRANSPOSE(KBLOCK)
!        PRINT *, KBEGNV
        IINVERSE(J1,1,1) = 1.0D0/KBEGNV(1)
        IINVERSE(J1,2,2) = 1.0D0/KBEGNV(2)
        IINVERSE(J1,3,3) = 1.0D0/KBEGNV(3)
        DO J2 = 1, NSITEPERBODY(J1)
           SITESRIGIDBODY(J2,:,J1) = MATMUL(KBLOCK,SITESRIGIDBODY(J2,:,J1))
        ENDDO
     ENDDO
  ENDIF
  
! make sure the two atoms used as reference for rigid bodies are suitable
! Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
  DO J1 = 1, NRIGIDBODY
     REFVECTOR(J1) = 1
     RTEST = .TRUE.
     DO WHILE (RTEST)
        RTEST = .FALSE.
        DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
           PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
           IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
              RTEST = .TRUE.
           ENDIF
        ENDDO
        PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
        IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
        IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
           RTEST = .TRUE.
        ENDIF
        IF (RTEST) THEN
           REFVECTOR(J1) = REFVECTOR(J1) + 1               
        ENDIF
     ENDDO
  ENDDO
END SUBROUTINE

!-----------------------------------------------------------
SUBROUTINE GENRIGID_READ_FROM_FILE ()
  USE INPUTmod, ONLY: L2U    
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  CHARACTER(LEN=10) CHECK1
  INTEGER :: J1, J2, DUMMY, iostatus
  REAL(KIND = REAL64) :: INICOORDS(3*NATOMS)

! read atomistic coordinates
! in future, no need for separate coordsinirigid
! currently the input coords files vary for CHARMM, AMBER, and RIGID BODIES
  IF (NATOMS == 0) THEN
     PRINT *, "ERROR STOP NOW > During generalised rigid body initialisation NATOMS = 0"
     STOP
  ENDIF

  OPEN(UNIT = 28, FILE = 'coordsinirigid', STATUS = 'OLD')
  DO J1 = 1, NATOMS
     READ(28, *) INICOORDS(3*J1-2), INICOORDS(3*J1-1), INICOORDS(3*J1)
  ENDDO
  CLOSE(UNIT = 28)

! determine no of rigid bodies
  NRIGIDBODY=0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO
     READ(222,*,IOSTAT=iostatus) CHECK1
     CALL L2U(CHECK1)
     IF (iostatus<0) THEN
        CLOSE(222)
        EXIT
     ELSE IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') then
        NRIGIDBODY=NRIGIDBODY+1
     ENDIF
  END DO
  CLOSE(222)
  
! determine maximum no of rigid body sites
  MAXSITE = 0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1,DUMMY
     IF (MAXSITE < DUMMY) MAXSITE = DUMMY
     DO J2 = 1, DUMMY
        READ(222,*) CHECK1
     ENDDO
  ENDDO
  CLOSE(222)

! Calling function for allocation to make more general (setup rigid bodies from code)
  CALL GENRIGID_ALLOCATE(NRIGIDBODY,MAXSITE)

! initialise SITESRIGIDBODY, RIGIDGROUPS, RIGIDSINGLES 
  OPEN(UNIT=222,FILE='rbodyconfig',status='unknown')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1, NSITEPERBODY(J1)
     DO J2 = 1, NSITEPERBODY(J1)
        READ(222,*) RIGIDGROUPS(J2,J1)
        ! check to make sure the current atom is not already in a rigid body
        IF (RIGIDISRIGID(RIGIDGROUPS(J2,J1))) THEN
            PRINT *," genrigid> ERROR: atom ",RIGIDGROUPS(J2,J1)," is in multiple rigid bodies! Stopping."
            STOP
        ELSE       
            ! if not, flag the current atom
            RIGIDISRIGID(RIGIDGROUPS(J2,J1))=.TRUE.
            RB_BY_ATOM(RIGIDGROUPS(J2,J1)) = J1
        ENDIF
        ! Moved initialization of coordinates to GENRIGID_INITIALISE, here only read the setup
        ! SITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY,1)
     ENDDO
  ENDDO
  CLOSE(222)
  CALL GENRIGID_INITIALISE(INICOORDS)
END SUBROUTINE GENRIGID_READ_FROM_FILE

!-----------------------------------------------------------

SUBROUTINE DEALLOCATE_GENRIGID()

  DEALLOCATE(NSITEPERBODY)
  DEALLOCATE(SITESRIGIDBODY)
  DEALLOCATE(RIGIDGROUPS)
  DEALLOCATE(REFVECTOR)
  DEALLOCATE(GR_WEIGHTS)
  DEALLOCATE(IINVERSE)
  DEALLOCATE(RIGIDISRIGID)
  DEALLOCATE(RB_BY_ATOM)
  DEALLOCATE(RIGIDSINGLES)
  DEALLOCATE(RIGIDCOORDS)

END SUBROUTINE DEALLOCATE_GENRIGID

!-----------------------------------------------------------

SUBROUTINE TRANSFORMRIGIDTOC (CMIN, CMAX, XCOORDS, XRIGIDCOORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J5, J7, J9
  INTEGER :: CMIN, CMAX  ! These will nearly always be 1,NRIGIDBODY respectively.
  REAL(KIND = REAL64) :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  REAL(KIND = REAL64) :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS)
  REAL(KIND = REAL64) :: COM(3) ! center of mass
  LOGICAL          :: GTEST !, ATOMTEST
  REAL(KIND = REAL64) :: MLATTICE(3,3)
  
  GTEST = .FALSE.

! no lattice - use identity matrix
  MLATTICE = 0D0
  MLATTICE(1,1)=1d0
  MLATTICE(2,2)=1D0
  MLATTICE(3,3)=1D0



  ! coord transformations for rigid bodies CMIN to CMAX
  DO J1 = CMIN, CMAX
     J5   = 3*J1
     J7   = 3*NRIGIDBODY + J5
     P(:) = XRIGIDCOORDS(J7-2:J7)
     CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! MLATTICE can have lattice transformation or be identity matrix
     COM = matmul(MLATTICE, XRIGIDCOORDS(J5-2:J5))
     DO J2 = 1,  NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        XCOORDS(3*J9-2:3*J9) = COM + MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
     ENDDO   
  ENDDO
  
! now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        J9 = RIGIDSINGLES(J1)
        XCOORDS(3*J9-2:3*J9) = XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1)
     ENDDO
  ENDIF
      
END SUBROUTINE TRANSFORMRIGIDTOC

!----------------------------------------------------------

SUBROUTINE ROTATEINITIALREF ()
IMPLICIT NONE
REAL(KIND = REAL64) :: P(3)
INTEGER J1

! rotate the system - new
  P(:) = OPTIMROTAVALUES(:)
!  P(1) = -1.0D0 * 8.0D0 * ATAN(1.0D0) 
!  P(2) = 0.0D0 !4.0D0 * ATAN(1.0D0) !-(8*ATAN(1.0D0) - 5.0D0)/DSQRT(2.0D0)
!  P(3) = 0.0D0 !4.0D0 * ATAN(1.0D0)
  DO J1 = 1, NRIGIDBODY
     CALL REDEFINERIGIDREF (J1,P)
  ENDDO

END SUBROUTINE ROTATEINITIALREF

!----------------------------------------------------------

SUBROUTINE REDEFINERIGIDREF (J1,P)

  IMPLICIT NONE
  
  INTEGER :: J1, J2     !No of processor
  REAL(KIND = REAL64) :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

!  PRINT *, "REDEFINE ", J1
  CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)  
  DO J2 = 1, NSITEPERBODY(J1)
     SITESRIGIDBODY(J2,:,J1) = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
  ENDDO

END SUBROUTINE REDEFINERIGIDREF

!----------------------------------------------------------

SUBROUTINE TRANSFORMCTORIGID (XCOORDS, XRIGIDCOORDS)
  USE COMMONS, ONLY: NATOMS, MYUNIT
  USE VEC3
  USE ROTATIONS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J9     !No of processor
  REAL(KIND = REAL64) :: COM(3), MASS
  REAL(KIND = REAL64) :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS)

! lattice matrix and inverse
  REAL(KIND = REAL64) MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

! extra variables for minpermdist
  REAL(KIND = REAL64) :: D, RMAT(3,3) 
  REAL(KIND = REAL64) :: PP1(3*NATOMS), PP2(3*NATOMS)
  
  XRIGIDCOORDS(:)=0.0D0
  NLATTICECOORDS=0
  MLATTICE=0
  MLATTICE(1,1)=1
  MLATTICE(2,2)=1
  MLATTICE(3,3)=1

  CALL INVERT3X3(MLATTICE, MLATTICEINV)

! loop over all rigid bodies
  DO J1 = 1, NRIGIDBODY
     COM = 0.0D0
     MASS = 0.0D0
     ! calculate center of mass
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        COM = COM + XCOORDS(3*J9-2:3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     COM = COM / MASS
     XRIGIDCOORDS(3*J1-2:3*J1) = COM

     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        PP1(3*J2-2:3*J2) = XCOORDS(3*J9-2:3*J9) - COM
        PP2(3*J2-2:3*J2) = SITESRIGIDBODY(J2,:,J1)
     ENDDO

     ! This could possibly be replaced with a call to ALIGN_DECIDE instead (might be faster and possibly slightly more
     ! accurate). Need to test this!
     CALL MINPERMDIST(PP1(1:3*NSITEPERBODY(J1)),PP2(1:3*NSITEPERBODY(J1)),NSITEPERBODY(J1),.FALSE.,D,RMAT)

     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = rot_mx2aa(RMAT)

     IF ( D/NSITEPERBODY(J1) > 0.1D0 ) THEN
        !print *, 'not going so well...'
        WRITE(MYUNIT, '(A, I3)')  'Warning: Genrigid > mapping looks bad for RB no ', J1 
        WRITE(MYUNIT, '(A)')  'Warning: Genrigid >  Often it is the permutation of the RB members, e.g. Hs in NH3'
     ENDIF

  ENDDO

! now translate everything to reduced units
  DO J1 = 1, NRIGIDBODY
    XRIGIDCOORDS(3*J1-2:3*J1) = MATMUL(MLATTICEINV, XRIGIDCOORDS(3*J1-2:3*J1))
  END DO

! now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        !added lattice stuff
        XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = MATMUL(MLATTICEINV, XCOORDS(3*J9-2:3*J9))
     ENDDO
  ENDIF

END SUBROUTINE TRANSFORMCTORIGID

!-----------------------------------------------------------

SUBROUTINE TRANSFORMGRAD (G, XR, GR)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  REAL(KIND = REAL64) :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  REAL(KIND = REAL64) :: PI(3)
  REAL(KIND = REAL64) :: DR1(3),DR2(3),DR3(3) 
  REAL(KIND = REAL64) :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  LOGICAL :: GTEST
  INTEGER :: NLATTICECOORDS
  
  NLATTICECOORDS=0

  GTEST = .TRUE.
  GR(:) = 0.0D0
  
  DO J1 = 1, NRIGIDBODY
     
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)

     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)

! translation
        GR(3*J1-2:3*J1) = GR(3*J1-2:3*J1) + G(3*J9-2:3*J9)
        
! rotation
        DR1(:) = MATMUL(DRMI1,SITESRIGIDBODY(J2,:,J1))
        DR2(:) = MATMUL(DRMI2,SITESRIGIDBODY(J2,:,J1))
        DR3(:) = MATMUL(DRMI3,SITESRIGIDBODY(J2,:,J1))
        GR(3*NRIGIDBODY+3*J1-2) = GR(3*NRIGIDBODY+3*J1-2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        GR(3*NRIGIDBODY+3*J1-1) = GR(3*NRIGIDBODY+3*J1-1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        GR(3*NRIGIDBODY+3*J1)   = GR(3*NRIGIDBODY+3*J1)   + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
  ENDDO


  IF (FREEZERIGIDBODYT .EQV. .TRUE.) THEN
     DO J1=1,NRIGIDBODY
        IF (FROZENRIGIDBODY(J1)) THEN
!           GR(3*J1-2:3*J1) = 0.0D0
!           GR(DEGFREEDOMS/2+3*J1-2:DEGFREEDOMS/2+3*J1) = 0.0D0
            CYCLE
        ENDIF
     ENDDO
!     GR(3*NRIGIDBODY-2:3*NRIGIDBODY) = 0.0D0
!     GR(6*NRIGIDBODY-2:6*NRIGIDBODY) = 0.0D0
  ENDIF

! single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY - NLATTICECOORDS) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = G(3*J9-2:3*J9)
     ENDDO
  ENDIF


END SUBROUTINE TRANSFORMGRAD

!-----------------------------------------------------------

SUBROUTINE AACONVERGENCE (G, XR, GR, RMS)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  REAL(KIND = REAL64) :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  REAL(KIND = REAL64) :: PI(3)
  REAL(KIND = REAL64) :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  REAL(KIND = REAL64) :: TORQUE(3)
  REAL(KIND = REAL64) :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  REAL(KIND = REAL64) :: DR1(3),DR2(3),DR3(3)
  real(kind = real64), intent(out) :: RMS 

  RMS = 0.0D0
  PI = (/0.0D0, 0.0D0, 0.0D0/)
  CALL RMDRVT(PI, RMI0, DRMI10, DRMI20, DRMI30, .TRUE.)

  DO J1 = 1, NRIGIDBODY
     
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)

     TORQUE(:) = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        DR1(:) = MATMUL(DRMI10,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR2(:) = MATMUL(DRMI20,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR3(:) = MATMUL(DRMI30,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        TORQUE(1) = TORQUE(1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        TORQUE(2) = TORQUE(2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        TORQUE(3) = TORQUE(3) + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
     TORQUE = MATMUL(TRANSPOSE(RMI), TORQUE)
     RMS = RMS + DOT_PRODUCT(TORQUE, MATMUL(TRANSPOSE(IINVERSE(J1,:,:)),TORQUE))
     RMS = RMS + 1.0D0/NSITEPERBODY(J1) * DOT_PRODUCT(GR(3*J1-2:3*J1),GR(3*J1-2:3*J1)) 
  ENDDO

  RMS=MAX(DSQRT(RMS/DEGFREEDOMS),1.0D-100)

END SUBROUTINE AACONVERGENCE

!--------------------------------------------------------------

!  Often we want to check if the atoms grouped in a rigid body has moved or not
!  They should not if everything is done correctly
!  REDEFINESITEST = .FALSE. then it prints to standard output
!  REDEFINESITEST = .TRUE. then regroup atoms, SITESRIGIDBODY rewritten

SUBROUTINE CHECKSITES (REDEFINESITEST, COORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  INTEGER :: J1, J2, DUMMY
  REAL(KIND = REAL64) :: XMASS, YMASS, ZMASS, PNORM, MASS
  REAL(KIND = REAL64) :: XSITESRIGIDBODY(MAXSITE,3,NRIGIDBODY)
  REAL(KIND = REAL64) :: COORDS(3*NATOMS)
  LOGICAL :: RTEST, REDEFINESITEST
  

  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY = RIGIDGROUPS(J2,J1)
        XSITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + XSITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XSITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XSITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        XSITESRIGIDBODY(J2,1,J1) = XSITESRIGIDBODY(J2,1,J1) - XMASS
        XSITESRIGIDBODY(J2,2,J1) = XSITESRIGIDBODY(J2,2,J1) - YMASS
        XSITESRIGIDBODY(J2,3,J1) = XSITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO
  ENDDO
  

  IF (REDEFINESITEST) THEN
!     PRINT *, " SITES REDEFINED "
     SITESRIGIDBODY(:,:,:) = XSITESRIGIDBODY(:,:,:)

!Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
     DO J1 = 1, NRIGIDBODY
        REFVECTOR(J1) = 1
        RTEST = .TRUE.
        DO WHILE (RTEST)
           RTEST = .FALSE.
           DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
              PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
              IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
                 RTEST = .TRUE.
              ENDIF
           ENDDO
           PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
           IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
           IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
              RTEST = .TRUE.
           ENDIF
           IF (RTEST) THEN
              REFVECTOR(J1) = REFVECTOR(J1) + 1               
           ENDIF
        ENDDO
     ENDDO
  ELSE
!     PRINT *, XSITESRIGIDBODY
  ENDIF

END SUBROUTINE CHECKSITES

!--------------------------------------------------------------

SUBROUTINE RBDET(A, DET)
  
  IMPLICIT NONE
  REAL(KIND = REAL64) :: A (3,3), DET

  DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)) 

END SUBROUTINE RBDET

SUBROUTINE INVERSEMATRIX(A, AINVERSE)
  
  IMPLICIT NONE
  REAL(KIND = REAL64) :: A (3,3), AINVERSE(3,3), DET

  DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1)) 
  AINVERSE(1,1) = A(2,2)*A(3,3)-A(2,3)*A(3,2)
  AINVERSE(1,2) = A(3,2)*A(1,3)-A(3,3)*A(1,2)
  AINVERSE(1,3) = A(2,3)*A(1,2)-A(2,2)*A(1,3)
  AINVERSE(2,1) = A(3,1)*A(2,3)-A(3,3)*A(2,1)
  AINVERSE(2,2) = A(1,1)*A(3,3)-A(1,3)*A(3,1)
  AINVERSE(2,3) = A(2,1)*A(1,3)-A(2,3)*A(1,1)
  AINVERSE(3,1) = A(3,2)*A(2,1)-A(3,1)*A(2,2)
  AINVERSE(3,2) = A(3,1)*A(1,2)-A(3,2)*A(1,1)
  AINVERSE(3,3) = A(2,2)*A(1,1)-A(2,1)*A(1,2)
  AINVERSE(:,:) = AINVERSE(:,:)/DET

END SUBROUTINE INVERSEMATRIX


! random rotation move for rigid bodies
SUBROUTINE GENRIGID_ROTATE(XCOORDS, ROTATEFACTOR)

USE COMMONS, ONLY: NATOMS
IMPLICIT NONE

INTEGER :: J1, J2, J3
REAL(KIND = REAL64) :: COM(3), MASS, PI, TWOPI, DPRAND
REAL(KIND = REAL64) :: ROTATIONMATRIX(3,3), TOROTATE(3), ATOMROTATED(3)
REAL(KIND = REAL64) :: RANDOMPHI, RANDOMTHETA, RANDOMPSI, ST, CT, SPH, CPH, SPS, CPS
REAL(KIND = REAL64), INTENT(INOUT) :: XCOORDS(3*NATOMS)
REAL(KIND = REAL64), INTENT(IN) :: ROTATEFACTOR

ROTATIONMATRIX(:,:) = 0.0D0
TOROTATE(:) = 0.0D0
! Define some constants
PI=ATAN(1.0D0)*4
TWOPI=2.0D0*PI

! Loop over all rigid bodies
DO J1 = 1, NRIGIDBODY
   IF (.NOT.FROZENRIGIDBODY(J1)) THEN
      COM = 0.0D0
      MASS = 0.0D0

! For each rigid body, calculate center of mass
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         COM = COM + XCOORDS(3*J3-2:3*J3)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
         MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
      ENDDO
      COM = COM / MASS

! Move the rigid body centre of mass to the origin
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) - COM
      ENDDO

! Calculate a random rotation scaled by ROTATEFACTOR
      RANDOMPHI=(DPRAND()-0.5)*TWOPI*ROTATEFACTOR
      RANDOMTHETA=(DPRAND()-0.5)*PI*ROTATEFACTOR
      RANDOMPSI=(DPRAND()-0.5)*TWOPI*ROTATEFACTOR
      ST=SIN(RANDOMTHETA)
      CT=COS(RANDOMTHETA)
      SPH=SIN(RANDOMPHI)
      CPH=COS(RANDOMPHI)
      SPS=SIN(RANDOMPSI)
      CPS=COS(RANDOMPSI)

! Assemble the rotation matrix
      ROTATIONMATRIX(1,1)=CPS*CPH-CT*SPH*SPS
      ROTATIONMATRIX(2,1)=CPS*SPH+CT*CPH*SPS
      ROTATIONMATRIX(3,1)=SPS*ST
      ROTATIONMATRIX(1,2)=-SPS*CPH-CT*SPH*CPS
      ROTATIONMATRIX(2,2)=-SPS*SPH+CT*CPH*CPS
      ROTATIONMATRIX(3,2)=CPS*ST
      ROTATIONMATRIX(1,3)=ST*SPH
      ROTATIONMATRIX(2,3)=-ST*CPH
      ROTATIONMATRIX(3,3)=CT

! Apply the rotation matrix to the atoms in the rigid body
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         TOROTATE = XCOORDS(3*J3-2:3*J3)
         ATOMROTATED=MATMUL(ROTATIONMATRIX,TOROTATE)
         XCOORDS(3*J3-2:3*J3) = ATOMROTATED
      ENDDO

! Translate the rigid body centre of mass back to its old position
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) + COM
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_ROTATE

! random rotation move for rigid bodies
SUBROUTINE GENRIGID_TRANSLATE(XCOORDS, TRANSLATEFACTOR)

USE COMMONS, ONLY: NATOMS
IMPLICIT NONE

INTEGER :: J1, J2, J3  
REAL(KIND = REAL64) DPRAND
REAL(KIND = REAL64), INTENT(INOUT) :: XCOORDS(3*NATOMS)
REAL(KIND = REAL64), INTENT(IN) :: TRANSLATEFACTOR
REAL(KIND = REAL64):: TRANSLATEVECTOR(3)

! Loop over all rigid bodies
DO J1 = 1, NRIGIDBODY
   IF (.NOT.FROZENRIGIDBODY(J1)) THEN
      DO J2=1,3
         TRANSLATEVECTOR(J2)=2.0*(DPRAND()-0.5)*TRANSLATEFACTOR
      ENDDO  
! Move the rigid body
      DO J2 = 1, NSITEPERBODY(J1)
         J3 = RIGIDGROUPS(J2, J1)
         XCOORDS(3*J3-2:3*J3) = XCOORDS(3*J3-2:3*J3) + TRANSLATEVECTOR
      ENDDO
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_TRANSLATE

! subroutine to update the reference coordinates for the rigid bodies using 
! the NATOMS coordinates in XCOORDS. Note that the rigid body coordinates are relative
! to the COM of each rigid body. 
SUBROUTINE GENRIGID_UPDATE_REFERENCE(XCOORDS)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2, DUMMY
  REAL(KIND = REAL64), INTENT(IN) :: XCOORDS(3*NATOMS)
  REAL(KIND = REAL64) :: XMASS, YMASS, ZMASS, MASS

  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0d0
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        XMASS = XMASS + XCOORDS(3*DUMMY-2)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XCOORDS(3*DUMMY-1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XCOORDS(3*DUMMY  )*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS  = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        SITESRIGIDBODY(J2,1,J1) = XCOORDS(3*DUMMY-2) - XMASS
        SITESRIGIDBODY(J2,2,J1) = XCOORDS(3*DUMMY-1) - YMASS
        SITESRIGIDBODY(J2,3,J1) = XCOORDS(3*DUMMY  ) - ZMASS
     ENDDO
  ENDDO

END SUBROUTINE GENRIGID_UPDATE_REFERENCE


! Subroutine to calculate the distance between the centre of mass of two rigid bodies
! XCOORDS must be in RIGID BODY FORMAT!
SUBROUTINE GENRIGID_COMDISTANCE(BODY1,BODY2,XCOORDS,COMDISTANCE)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: BODY1, BODY2
  REAL(KIND = REAL64), INTENT(IN) :: XCOORDS(3*NATOMS)
  REAL(KIND = REAL64), INTENT(OUT) :: COMDISTANCE
  REAL(KIND = REAL64) :: BODY1COM(3), BODY2COM(3)

! For each rigid body, calculate center of mass
! Body 1
BODY1COM(:) = 0.0D0
BODY2COM(:) = 0.0D0

! You MUST be using rigid body coordinates here
! Assign BODY1/BODY2 COMs
BODY1COM = XCOORDS(3*BODY1-2:3*BODY1)
BODY2COM = XCOORDS(3*BODY2-2:3*BODY2)

! Calculate the COM-COM distance
COMDISTANCE = SQRT((BODY1COM(1)-BODY2COM(1))**2+(BODY1COM(2)-BODY2COM(2))**2+(BODY1COM(3)-BODY2COM(3))**2)


END SUBROUTINE GENRIGID_COMDISTANCE

! Subroutine to check that all rigid bodies are within DIST of each other. Returns TEST=.TRUE. if they are.
! We only require each rigid body to be within DIST of one other body to allow for large homogeneous systems.
! This might need changing in the future to require a user defined number of bodies to be within DIST.
SUBROUTINE GENRIGID_DISTANCECHECK(XCOORDS,DIST,TEST)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2
  REAL(KIND = REAL64), INTENT(IN) :: XCOORDS(3*NATOMS), DIST
  REAL(KIND = REAL64) :: COMDIST, MINCOMDIST
  LOGICAL, INTENT(OUT) :: TEST

TEST = .TRUE.
! Loop over all pairs of rigid bodies
DO J1 = 1,NRIGIDBODY
! Set the initial minimum distance to a large value
   MINCOMDIST = HUGE(1.0D0)
   DO J2 = 1,NRIGIDBODY
! Skip self-self checks
      IF (J1.EQ.J2) CYCLE
! Check the COM-COM distance
      CALL GENRIGID_COMDISTANCE(J1,J2,XCOORDS,COMDIST)
! If the calculated distance is less than the current minimum, replace it
      IF (COMDIST.LT.MINCOMDIST) THEN
         MINCOMDIST = COMDIST
      ENDIF
   ENDDO
! Once body J1 has been compare to all bodies J2, check at least one is within DIST
   IF (MINCOMDIST.GT.DIST) THEN
      TEST = .FALSE.
      RETURN
   ENDIF
ENDDO

END SUBROUTINE GENRIGID_DISTANCECHECK

!----------------------------------------------------------------------------------------------
! See section 2.2.1 of Mochizuki et al, PCCP 16 (2014), 2842-53  DOI 10.1039/C3CP53537A
! However, the notation used is closer to Wales and Ohmine, JCP 98 (1993), 7257-7268 (KBLOCK, U)
! See also rigidb.f90>NRMLMD, an obsolete version of the same procedure which has comments from the
! original authors. The comments on this routine were mostly added by sn402.

SUBROUTINE GENRIGID_EIGENVALUES(X, ATOMMASS, DIAG, INFO)
  
  USE COMMONS
  USE MODHESS
  IMPLICIT NONE

  INTEGER          :: IR, IC, OFFSET, NDIM
  INTEGER          :: I, J, J1, J2, J3, J5, J8, K1, K2, ISTART, IFINISH, JSTART, JFINISH 
  REAL(KIND = REAL64) :: U(DEGFREEDOMS,DEGFREEDOMS), KBLOCK(3,3), KBEGNV(3), TMASS, ENERGY
  REAL(KIND = REAL64) :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS), G(3*NATOMS)
  REAL(KIND = REAL64) :: X(3*NATOMS), FRQN(DEGFREEDOMS), ATOMMASS(NATOMS), DIAG(3*NATOMS)
  REAL(KIND = REAL64) :: P(3), RMI(3,3), DRMI(3,3), DR(3)
  REAL(KIND = REAL64) :: KD(DEGFREEDOMS), AP(DEGFREEDOMS,DEGFREEDOMS)
! the following required to call the LAPACK routine DSYEV
  INTEGER          :: INFO
  INTEGER, PARAMETER :: LWORK = 1000000 ! the dimension is set arbitrarily
  REAL(KIND = REAL64) :: WORK(LWORK)
  LOGICAL          :: ART

  ! The METRICTENSOR keyword instructs us to override all calls to this subroutine with a call
  ! to the newer subroutine GENRIGID_NORMALMODES. They should have exactly the same behaviour.
! IF(METRICTENSOR) THEN
!    CALL GENRIGID_NORMALMODES(X, ATOMMASS, DIAG, INFO)
!    RETURN
! ENDIF

! Initialize

  OFFSET = 3*NRIGIDBODY
  U(:,:) = 0.D0
  IR     = 0  ! IR+1 is the index of the first coordinate for the current rigid body
  IC     = 0
  ART = ATOMRIGIDCOORDT

! Transform coordinates

  IF ( ATOMRIGIDCOORDT .EQV. .TRUE. ) THEN
     CALL TRANSFORMCTORIGID (X, XRIGIDCOORDS)
  ELSE
     XRIGIDCOORDS(1:DEGFREEDOMS) = X(1:DEGFREEDOMS)
  ENDIF

! OPEN(UNIT=1356,FILE='rigidcoords',STATUS='OLD')
! PRINT *,'rigid coords, NRIGIDBODY=',NRIGIDBODY
! DO J1=1,6*NRIGIDBODY
! READ(1356,*) XRIGIDCOORDS(J1)
! PRINT *,J1,XRIGIDCOORDS(J1)
! ENDDO
! CLOSE(1356)

  DO J1 = 1, NRIGIDBODY

     J3 = 3*J1
     J5 = OFFSET + J3
     P  = XRIGIDCOORDS(J5-2:J5)
     ! KBLOCK is the moment of inertia matrix for this rigid body.
     KBLOCK(:,:) = 0.D0     
     ! Get the rotation matrix RMI that corresponds to P (doesn't calculate any derivatives)
     CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)

     TMASS = 0.0D0 ! The total mass of the rigid body
     DO J2 = 1, NSITEPERBODY(J1)
        J8 = RIGIDGROUPS(J2, J1)
        ! DR is the coordinates of this atom relative to the rigid body CoM.
        DR(:)  = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
        TMASS = TMASS + ATOMMASS(J8)
        DO I = 1, 3
           ! Diagonal terms in the moment of inertia about the rigid-body CoM
           ! sn402: Is this right? Shouldn't Ixx = sum_i m_i * (y_i^2+z_i^2) ? Here we have Ixx = sum_i m_i * (x_i^2+y_i^2+z_i^2)
           KBLOCK(I,I) = KBLOCK(I,I) + ATOMMASS(J8)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
           DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
              KBLOCK(I,J) = KBLOCK(I,J) - ATOMMASS(J8)*DR(I)*DR(J)
           ENDDO
        ENDDO
     ENDDO
     ! Diagonalise KBLOCK
     ! KBEGNV are the KBLOCK eigenvalues, KBLOCK now contains the eigenvectors
     CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
     IF (INFO /= 0) THEN
        WRITE(*,*) 'NRMLMD > Error in DSYEV with KBLOCK, INFO =', INFO
        STOP
     ENDIF

!    Construction of the matrix U, which diagonalises the coordinate vector x = (r,p) for this rigid body.
!    The reason we want this is that making the transform w = matmul(U,p) gives us coordinates
!    which are diagonal in the kinetic energy.

!    First: translation coordinates (recall IR+1 is the first coordinate for this rigid body)
!    No diagonalisation required. U is identity.
     U(IR+1,IC+1) = 1.D0; U(IR+2,IC+2) = 1.D0; U(IR+3,IC+3) = 1.D0
     ! KD is the diagonalised KBLOCK. Only diagonal elements are recorded.
     ! Scaling by sqrt(TMASS) is required to move to the CoM frame, as usual in a normal mode calculation.
     KD(IC+1:IC+3) = 1.D0/SQRT(TMASS)

!    Now rotational coordinates.
!    Recall that KBLOCK now contains the eigenvectors of the original matrix. These are required to diagonalise
!    the rotational component of the coordinate vector.
     U(OFFSET+IR+1:OFFSET+IR+3,OFFSET+IC+1:OFFSET+IC+3) = KBLOCK(:,:)
     KD(OFFSET+IC+1:OFFSET+IC+3) = 1.D0/SQRT(KBEGNV(:))               ! See Mochizuki et al

     ! IR = IC = 3*(J1-1) - these two variables indicate the first coordinate position of the current rigid body.
     IR = IR + 3
     IC = IC + 3 

  ENDDO ! Loop over rigid bodies (J1)

  ! Deal with the free atoms. No off-diagonal terms in the kinetic energy here, so they are automatically diagonalised and only
  ! require mass-rescaling.
  DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J8 = RIGIDSINGLES(J1)
     U(6*NRIGIDBODY+3*J1-2,6*NRIGIDBODY+3*J1-2) = 1.D0
     U(6*NRIGIDBODY+3*J1-1,6*NRIGIDBODY+3*J1-1) = 1.D0
     U(6*NRIGIDBODY+3*J1  ,6*NRIGIDBODY+3*J1  ) = 1.D0
     KD(6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1 ) = 1.D0/SQRT(ATOMMASS(J8))     
  ENDDO

  ! We have now obtained the required diagonal coordinates for the total kinetic energy.
  ! We next need to obtain the Hessian in terms of these coordinates.

  RBAANORMALMODET = .TRUE.
  ATOMRIGIDCOORDT = .FALSE.
  XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
  XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
  ! STEST is set to TRUE so we get the Hessian calculated.
  ! When TRANSFORM_HESSIAN is called from POTENTIAL, it will behave differently because of RBAANORMALMODET=TRUE.
  ! See TRANSFORM_HESSIAN for comments.
  CALL POTENTIAL(XCOORDS,G,ENERGY,.TRUE.,.TRUE.)

  RBAANORMALMODET = .FALSE.
!  ATOMRIGIDCOORDT = .TRUE.

  NDIM = DEGFREEDOMS  
  AP(:,:) = 0.D0  ! This will be our Hessian in diagonalised general coordinates.
  ! Fill in the lower-left triangle of the transformed Hessian
  DO I = 1, NDIM
     DO J = 1, I
        
        ! sn402: I'm not quite sure what's going on here.
        ! This block seems to suggest that the block AP(1:3,1:3) gets written to more than everything else.
        ! Wouldn't it be easier to use the upper-right triangle and use
        !DO I=1,NDIM
        !   DO J=I,NDIM
        !      DO K1=1,NDIM
        !         DO K2=1,NDIM
        !            ?
        ! Possibly this way is more efficient because U is 0 away from the main diagonal?
        IF (I .LE. 2) THEN
           ISTART = 1
        ELSE
           ISTART = I-2
        ENDIF
        IF (I .GE. NDIM -2) THEN
           IFINISH = NDIM
        ELSE
           IFINISH = I+2
        ENDIF
        IF (J .LE. 2) THEN
           JSTART = 1
        ELSE
           JSTART = J-2
        ENDIF
        IF (J .GE. NDIM-2) THEN
           JFINISH = NDIM
        ELSE
           JFINISH = J+2
        ENDIF

        DO K1 = ISTART, IFINISH
           DO K2 = JSTART, JFINISH
              ! The actual coordinate transform. Recall that U is diagonal for translational coordinates and single
              ! atoms (the elements are square roots of the corresponding masses)
              ! But U and HESS are both different for the rotational coordinates.
              AP(I,J) = AP(I,J) + U(K1,I)*HESS(K1,K2)*U(K2,J)
           ENDDO
        ENDDO
        AP(I,J) = KD(I)*AP(I,J)*KD(J) ! Weight the coordinates by the eigenvalues of the inertia tensor (see above)
     ENDDO
  ENDDO

  ! Now diagonalise the modified Hessian. The eigenvalues (returned as FRQN) are the normal mode frequencies.
  ! If DUMPV is set, then on return AP contains the eigenvectors of the transformed Hessian - i.e. the normal modes.
! IF (DUMPV) THEN
!    CALL DSYEV('V','L',NDIM,AP(1:NDIM,1:NDIM),NDIM,FRQN,WORK,LWORK,INFO)
! ELSE
     CALL DSYEV('N','L',NDIM,AP(1:NDIM,1:NDIM),NDIM,FRQN,WORK,LWORK,INFO)
! ENDIF

!  call eigensort_val_asc(FRQN,AP,NDIM,NDIM)
! Mass-weight the normal mode frequencies here?
! Looks like someone has hard-coded the units in. Need to change that.
  IF ((.TRUE.)) THEN ! .OR. CHARMMT)) THEN  ! I need to get this included at some point
     ! Hessian diagonalisation returns square frequencies in internal units. For AMBER/CHARMM, these units are
     ! (kCal mol^-1)/(amu Angstrom^2) = 4.184E26 s^-2.
     ! So to convert to rad s^-1, we want sqrt(FRQN(I)*4.184D26) = 2.045E13*sqrt(FRQN(I))
     ! For Hz: sqrt(FRQN(I)*4.184D26)/2*pi = 3.255E12*sqrt(FRQN(I))
     ! For cm^-1: sqrt(FRQN(I)*4.184D26)/(2*pi*c*100) = 108.52*sqrt(FRQN(I))
     DO I = 1, NDIM
        IF (FRQN(I) > 0.0D0) THEN
           FRQN(I) = SQRT((FRQN(I)))*108.52
        ELSE
           FRQN(I) = -SQRT((-FRQN(I)))*108.52
        ENDIF
     ENDDO
  ENDIF

  DIAG(1:DEGFREEDOMS) = FRQN
  DIAG(DEGFREEDOMS+1:3*NATOMS) = 0.0D0

! set the Hessian matrix to its eigenvector matrix
  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = AP
! Restore this variable to its saved value.
  ATOMRIGIDCOORDT = ART

!  OPEN(UNIT = 28, FILE = 'LRBNORMALMODES')
!  WRITE(28, *) ENERGY, 3*NATOMS, DEGFREEDOMS
!  DO J1 = 3, DEGFREEDOMS/3
!     WRITE(28, *) FRQN(3*J1-2), FRQN(3*J1-1), FRQN(3*J1)
!  ENDDO
!  CLOSE(UNIT = 28)

! IF (DUMPV) THEN      
!    CALL GENRIGID_EIGENVECTORS(AP, XRIGIDCOORDS, ATOMMASS)
! ENDIF
  
END SUBROUTINE GENRIGID_EIGENVALUES

!----------------------------------------------------------------------------------------------

SUBROUTINE GENRIGID_EIGENVECTORS(AP, XRIGIDCOORDS, ATOMMASS)

  USE COMMONS
  USE MODHESS
  IMPLICIT NONE
        
  INTEGER :: OFFSET, SMODE, I, J, J1, J2, J3, J5, J8
  REAL(KIND = REAL64) :: XTEMP(DEGFREEDOMS), XXTEMP(DEGFREEDOMS), AP(DEGFREEDOMS,DEGFREEDOMS)
  REAL(KIND = REAL64) :: X(3*NATOMS), XCOORDS(3*NATOMS), XRIGIDCOORDS(3*NATOMS)
  REAL(KIND = REAL64) :: P(3), KBLOCK(3,3), KBEGNV(3), RMI(3,3), RRMI(3,3), DRMI(3,3), DR(3)
  REAL(KIND = REAL64) :: TMASS, ATOMMASS(3*NATOMS)
  INTEGER, PARAMETER :: LWORK = 10000 ! the dimension is set arbitrarily
  REAL(KIND = REAL64) :: WORK(LWORK)
  INTEGER          :: INFO

  OFFSET = 3*NRIGIDBODY
  CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, X, XRIGIDCOORDS)

! Here are the normal modes you want to be computed - the first six are zero
  DO SMODE = 1, DEGFREEDOMS
     DO I = 1, DEGFREEDOMS
! This gives you DX where DX is taken along the direction of an eigenvector
        XTEMP(I) = AP(I,SMODE) * 0.05
     ENDDO
     DO J1 = 1, NRIGIDBODY
        J3 = 3*J1
! The rotational part - undiagonalized it
        J5 = OFFSET + J3
        P  = XRIGIDCOORDS(J5-2:J5)
        KBLOCK(:,:) = 0.D0                 
! Computes the rotation matrix which brings the reference geometry in stationary frame 
! to the atom positions in current minimum geometry
        CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)
                 
! Computing inertia matrix in the moving frame, i.e. current minimum geometry
        TMASS = 0.0D0
        DO J2 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J2, J1)        
           DR(:)  = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
           TMASS = TMASS + ATOMMASS(J8)
           DO I = 1, 3
              KBLOCK(I,I) = KBLOCK(I,I) + ATOMMASS(J8)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - ATOMMASS(J8)*DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
! Diagonalise inertia matrix
        CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)

        XTEMP(J5-2) = XTEMP(J5-2)/SQRT(KBEGNV(1)) 
        XTEMP(J5-1) = XTEMP(J5-1)/SQRT(KBEGNV(2)) 
        XTEMP(J5  ) = XTEMP(J5  )/SQRT(KBEGNV(3)) 
                 
! Going from the diagonalised rotation coordinates to per rigid body angle-axis coordinates in the moving frame
        XXTEMP(J5-2) = KBLOCK(1,1)*XTEMP(J5-2) + KBLOCK(1,2)*XTEMP(J5-1) + KBLOCK(1,3)*XTEMP(J5  )
        XXTEMP(J5-1) = KBLOCK(2,1)*XTEMP(J5-2) + KBLOCK(2,2)*XTEMP(J5-1) + KBLOCK(2,3)*XTEMP(J5  )
        XXTEMP(J5  ) = KBLOCK(3,1)*XTEMP(J5-2) + KBLOCK(3,2)*XTEMP(J5-1) + KBLOCK(3,3)*XTEMP(J5  )
                 
! Computing the rotation matrix which rotates the current minimum geometry due to its normal mode displacements
        P = XXTEMP(J5-2:J5)
        CALL RMDFAS(P, RRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)

! Compute the displaced positions of the atoms
        DO J2 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J2, J1)        
           XCOORDS(3*J8-2:3*J8) = XRIGIDCOORDS(J3-2:J3) + XTEMP(J3-2:J3)/SQRT(TMASS) &
                + MATMUL(RRMI(:,:),MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1)))
        ENDDO
     ENDDO

     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        J8 = RIGIDSINGLES(J1)
        J5 = 2*OFFSET + 3*J1
        XCOORDS(3*J8-2:3*J8) = XRIGIDCOORDS(J5-2:J5) + XTEMP(J5-2:J5)/SQRT(ATOMMASS(J8))
     ENDDO

! Computes the atoms eigenvectors
     DO I = 1, 3*NATOMS
        HESS(I,SMODE) = (XCOORDS(I)-X(I)) / 0.05
     ENDDO

  ENDDO

END SUBROUTINE GENRIGID_EIGENVECTORS

!-----------------------------------------------------------
! Calculate the Hessian in rigid body coordinates from the Gradient and Hessian in Cartesians.
! Follows the procedure outlined in Mochizuki et al, PCCP 16 (2014)
SUBROUTINE TRANSFORMHESSIAN (H, G, XR, HR, RBAANORMALMODET)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J3, J4, J8, J9, K, L
  REAL(KIND = REAL64), INTENT(IN) :: G(3*NATOMS), H(3*NATOMS,3*NATOMS), XR(DEGFREEDOMS)
  REAL(KIND = REAL64), INTENT(OUT) :: HR(DEGFREEDOMS,DEGFREEDOMS)
  REAL(KIND = REAL64) :: PI(3)
  REAL(KIND = REAL64) :: AD2R11(3),AD2R22(3),AD2R33(3),AD2R12(3),AD2R23(3),AD2R31(3) 
  REAL(KIND = REAL64) :: ADR1(3),ADR2(3),ADR3(3) 
  REAL(KIND = REAL64) :: ARMI(3,3), ADRMI1(3,3), ADRMI2(3,3), ADRMI3(3,3)
  REAL(KIND = REAL64) :: AD2RMI11(3,3), AD2RMI22(3,3), AD2RMI33(3,3)
  REAL(KIND = REAL64) :: AD2RMI12(3,3), AD2RMI23(3,3), AD2RMI31(3,3)
  REAL(KIND = REAL64) :: BDR1(3),BDR2(3),BDR3(3) 
  REAL(KIND = REAL64) :: BRMI(3,3), BDRMI1(3,3), BDRMI2(3,3), BDRMI3(3,3)
  REAL(KIND = REAL64) :: BD2RMI11(3,3), BD2RMI22(3,3), BD2RMI33(3,3)
  REAL(KIND = REAL64) :: BD2RMI12(3,3), BD2RMI23(3,3), BD2RMI31(3,3)
  LOGICAL :: GTEST, STEST, RBAANORMALMODET
  REAL(KIND = REAL64) :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  REAL(KIND = REAL64) :: D2RMI10(3,3), D2RMI20(3,3), D2RMI30(3,3), D2RMI120(3,3), D2RMI230(3,3), D2RMI310(3,3)
  
  GTEST = .TRUE.
  STEST = .TRUE.
  HR(:,:) = 0.0D0


  IF ( RBAANORMALMODET ) THEN
     PI = (/0.0D0, 0.0D0, 0.0D0/)
     CALL RMDFAS(PI, RMI0, DRMI10, DRMI20, DRMI30, D2RMI10, D2RMI20, D2RMI30, D2RMI120, D2RMI230, D2RMI310, GTEST, STEST)
  ENDIF

  DO J1 = 1, NRIGIDBODY

     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)
    
     DO J2 = J1, NRIGIDBODY

        PI = XR(3*NRIGIDBODY+3*J2-2 : 3*NRIGIDBODY+3*J2)
        CALL RMDFAS(PI, BRMI, BDRMI1, BDRMI2, BDRMI3, BD2RMI11, BD2RMI22, BD2RMI33, BD2RMI12, BD2RMI23, BD2RMI31, GTEST, STEST)

             
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)

           DO J4 = 1, NSITEPERBODY(J2)
              J9 = RIGIDGROUPS(J4, J2)

!  translation
              HR(3*J1-2:3*J1, 3*J2-2:3*J2) = HR(3*J1-2:3*J1, 3*J2-2:3*J2) + H(3*J8-2:3*J8, 3*J9-2:3*J9)

!  rotations
              IF ( RBAANORMALMODET ) THEN
                 ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 BDR1(:) = MATMUL(DRMI10,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR2(:) = MATMUL(DRMI20,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR3(:) = MATMUL(DRMI30,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
              ELSE
                 ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
                 ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
                 ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
                 BDR1(:) = MATMUL(BDRMI1,SITESRIGIDBODY(J4,:,J2))
                 BDR2(:) = MATMUL(BDRMI2,SITESRIGIDBODY(J4,:,J2))
                 BDR3(:) = MATMUL(BDRMI3,SITESRIGIDBODY(J4,:,J2))
              ENDIF

!  mixed translation rotation
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*J1  , 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*J1  , 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*J1-2, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*J1-1, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*J1  , 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR3(:))        
              
              IF (J2 > J1) THEN
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-2) = HR(3*J2  , 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-1) = HR(3*J2  , 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1  ) = HR(3*J2  , 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
              ENDIF
!  double rotation
              DO K = 1, 3
                 DO L = 1, 3
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR3(L)   
                 ENDDO
              ENDDO
           ENDDO ! Loop over J4

           IF (J1 .EQ. J2) THEN
              IF ( RBAANORMALMODET ) THEN
                 AD2R11(:) = MATMUL(D2RMI10, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R22(:) = MATMUL(D2RMI20, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R33(:) = MATMUL(D2RMI30, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R12(:) = MATMUL(D2RMI120,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R23(:) = MATMUL(D2RMI230,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R31(:) = MATMUL(D2RMI310,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ELSE
                 AD2R11(:) = MATMUL(AD2RMI11,SITESRIGIDBODY(J3,:,J1))
                 AD2R22(:) = MATMUL(AD2RMI22,SITESRIGIDBODY(J3,:,J1))
                 AD2R33(:) = MATMUL(AD2RMI33,SITESRIGIDBODY(J3,:,J1))
                 AD2R12(:) = MATMUL(AD2RMI12,SITESRIGIDBODY(J3,:,J1))
                 AD2R23(:) = MATMUL(AD2RMI23,SITESRIGIDBODY(J3,:,J1))
                 AD2R31(:) = MATMUL(AD2RMI31,SITESRIGIDBODY(J3,:,J1))
              ENDIF

              ! p_x, p_x
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R11(:))
              ! p_y, p_y
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R22(:))
              ! p_z, p_z
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R33(:))
              ! p_x, p_y
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R12(:))
              ! p_y, p_z
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R23(:))
              ! p_z, p_x
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R31(:))
           ENDIF
        ENDDO ! Loop over J3 (i.e. J8)

        IF (J1 .EQ. J2) THEN
           HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1-1) 
           HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1  )
           HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-2)
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
        ELSE
           HR(3*J2-2:3*J2, 3*J1-2:3*J1) = TRANSPOSE(HR(3*J1-2:3*J1, 3*J2-2:3*J2))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J2-2:3*J2) = &
                TRANSPOSE(HR(3*J2-2:3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
                TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
        ENDIF

     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY

     DO J2 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
            
        J9 = RIGIDSINGLES(J2)
        PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
        CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)

        
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)
                     
!  translation
           HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) &
                + H(3*J8-2:3*J8, 3*J9-2:3*J9)

!  rotations
           IF ( RBAANORMALMODET ) THEN
              ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
           ELSE
              ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
              ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
              ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
           ENDIF
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
        ENDDO
        
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
             TRANSPOSE(HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
             TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
     ENDDO
  ENDDO

  DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J8 = RIGIDSINGLES(J1)
     DO J2 = J1, (DEGFREEDOMS - 6*NRIGIDBODY)/3            
        J9 = RIGIDSINGLES(J2)                     
        HR(6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = H(3*J8-2:3*J8, 3*J9-2:3*J9)
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1) = H(3*J9-2:3*J9, 3*J8-2:3*J8)
     ENDDO
  ENDDO

  !  safety check
! IF(DEBUG) THEN
!     DO J1 = 1,DEGFREEDOMS
!        DO J2 = 1,DEGFREEDOMS
!            IF(ABS(HR(J1,J2)-HR(J2,J1)) .GT. 1.0E-7) THEN
!                write(*,*) "transformHessian> Asymmetric Hessian, coords ", J1, J2
!            ENDIF
!        ENDDO
!     ENDDO
! ENDIF

END SUBROUTINE TRANSFORMHESSIAN

END MODULE GENRIGID

