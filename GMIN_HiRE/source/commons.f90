! Global variable definitions, module available to all parts of code

      MODULE COMMONS
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_PTR
      USE PREC, ONLY: REAL64
      IMPLICIT NONE
      SAVE

      TYPE NUCLEOTIDE
        INTEGER :: FATOM
        INTEGER :: LATOM
      END TYPE NUCLEOTIDE

      INTEGER :: MYNODE, MYUNIT, MYEUNIT, MYMUNIT, MYBUNIT !globals for output units and computing nodes
      INTEGER :: MYSEUNIT, MYSMUNIT                        !output for SAXS
      INTEGER :: NATOMS, NATOMSALLOC                       !number of atoms (need two variables one for use, one for allocation)
      INTEGER :: NPAR =1                                   !number of parallel runs 
      INTEGER :: NSAVE = 5                                 !number of saved structures
      INTEGER :: NACCEPT = 50                              !number of accepted steps
      INTEGER :: MCSTEPS = 10000                           !number of steps
      INTEGER :: MAXIT = 500 , MAXIT2 = 500                !maximum iteration for the sloppy and tight quenching
      INTEGER :: MUPDATE = 4                               !number of stored previous steps in LBFGS minimiser
      INTEGER :: NPCALL = 0                                !number of potential calls
      INTEGER :: CHECKDID = 1                              !ID for CHECKD keyword (0 is potential, 1 is gradient, 2 is Hessian, 3 is a single quench)
      INTEGER :: DUMPINT = 1000                            !interval for dumping 'GMIN.dump'
      INTEGER :: NSEED = 0                                 !random number seed?
      INTEGER :: QUENCHFRQ = 1                             !frequency of quenching (should be done every step?)
      INTEGER :: PRTFRQ = 1                                !print frequency for febh routines

      REAL(KIND = REAL64) :: TSTART                        !initial CPU time
      REAL(KIND = REAL64) :: CQMAX = 1.0D-8                !tight convergence limit
      REAL(KIND = REAL64) :: BQMAX = 1.0D-3                !sloppy convergence limit
      REAL(KIND = REAL64) :: RADIUS = 0.0D0                !radius of particle container
      REAL(KIND = REAL64) :: TFAC = 1.0D0                  !Temperature scaling factor	
      REAL(KIND = REAL64) :: GEOMDIFFTOL = 0.5D0           !Tolerance on distance for two structures to be classed as the same
      REAL(KIND = REAL64) :: SRATIO = -1.0D0               !Step size ratio for scaling
      REAL(KIND = REAL64) :: TRATIO = -1.0D0               !Temperature ratio for scaling
      REAL(KIND = REAL64) :: MAXERISE = 1.0D-10            !Maximum energy increase in minimisation
      REAL(KIND = REAL64) :: MAXBFGS = 0.4D0               !Maximum step in minimiser
      REAL(KIND = REAL64) :: IDIFFTOL = 0.1D0              !Tolerance on inertia (not used?)
      REAL(KIND = REAL64) :: ECONV = 0.02D0                !Tolerance on energy difference for two structures to be the same
      REAL(KIND = REAL64) :: DGUESS = 0.1D0                !Initial guess for Hessian elements
      REAL(KIND = REAL64) :: COLDFUSIONLIMIT = -1.0D6      !Energy limit considered as cold fusion
      REAL(KIND = REAL64) :: MAXEFALL = -HUGE(1.0D0)       !Maximum decrease in energy allowed
      REAL(KIND = REAL64) :: SYMFCTR = 0.0D0               ! 
      REAL(KIND = REAL64) :: GMAX  = 1.0D-8                !Current convergence (set to either BQMAX or CQMAX in quench.F)
      REAL(KIND = REAL64) :: POTEL = 0.0D0                 !Current potential energy (removed common block)
      REAL(KIND = REAL64) :: RMS = 0.0D0                   !Current RMS force
      REAL(KIND = REAL64) :: QTESTMAX = 1.0D-2             !convergence limit for CHECKDID 3


      LOGICAL :: DEBUG = .FALSE.                           !Enable debug printing and extra checks
      LOGICAL :: MPIT = .FALSE.                            !PT run using MPI
      LOGICAL :: RESTORET = .FALSE.                        !Restore previous run
      LOGICAL :: RANSEEDT = .FALSE.                        !Seed random numbers from input integer
      LOGICAL :: RATIOT = .FALSE.                          !Set ratio for temperature and step size
      LOGICAL :: MAXERISE_SET = .FALSE.                    !Set maximum energy increase
      LOGICAL :: TRACKDATAT = .FALSE.                      !Turn on tracking for additional properties  
      LOGICAL :: DUMPSTRUCTURES = .FALSE.                  !Dump pdb files
      LOGICAL :: DUMPMINT = .FALSE.                        !additional output
      LOGICAL :: SAVEQ = .TRUE.                            !save results
      LOGICAL :: USEFRQS = .FALSE.                         !Use quantum partition function
      LOGICAL :: COLDFUSION = .FALSE.                      !Has coldfusion been diagnosed?
      LOGICAL :: GRADPROBLEMT = .FALSE.                    !Any gradient probelsm encountered?
      LOGICAL :: CHECKDT  = .FALSE.                        !Run CHECKD routines
      LOGICAL :: DOCARTSTEP = .FALSE.                      !Do Cartesian step in this BH step?
 
      CHARACTER(LEN=130) :: DUMPFILE = ''

! Run-time variables, dimension will be NPAR to account for parallel runs
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NQ                             !Number of quenches for each run 

      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:,:) :: COORDS, COORDSO  !Coordinatees
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: ATMASS             !Mass of atoms/beads
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: TEMP               !Temperature
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: STEP               !Step size taken
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: OSTEP              !
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: ASTEP              !
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: ACCRAT             !
      REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: EPREV              !

      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FIXBOTH                        !
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FIXTEMP                        !
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FIXSTEP                        !
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: TMOVE                          !
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: OMOVE                          !

! Dumping Markov 
      LOGICAL :: DUMP_MARKOV=.FALSE.
      INTEGER :: DUMP_MARKOV_NWAIT=0
      INTEGER :: DUMP_MARKOV_NFREQ=1
      INTEGER, ALLOCATABLE :: DUMPXYZUNIT(:)

! Ratio variables
      REAL(KIND = REAL64) :: SUMSTEP = 0.D0, SUMTEMP = 0.D0

! Acceptance ratio method for step-size scaling
      LOGICAL :: ARMT = .FALSE.
      REAL(KIND = REAL64) :: ARMA = 0.4, ARMB = 0.4

! Q module variables necessary
      INTEGER, ALLOCATABLE :: FF(:)                                        !Number of steps when first found
      INTEGER, ALLOCATABLE :: NPCALL_QMIN(:)                               !Number of potential calls when found
      REAL(KIND = REAL64), ALLOCATABLE :: QMIN(:)                          !Energy of saved minima
      REAL(KIND = REAL64), ALLOCATABLE :: QMINP(:,:)                       !coordinates of saved minima
      INTEGER, ALLOCATABLE :: QMINT(:,:), QMINNATOMS(:)


! PT basin-hopping
      REAL(KIND = REAL64) ::  PTTMIN = 0.0D0, PTTMAX = 1.0D0, EXCHPROB = 0.2D0
      INTEGER :: EXCHINT = 100
      LOGICAL :: PTRANDOM=.FALSE., PTINTERVAL=.FALSE., PTSINGLE=.FALSE., PTSETS=.FALSE.

! set centre
      LOGICAL :: CENT = .FALSE., QUCENTRET = .FALSE., SETCENT = .FALSE.
      REAL(KIND = REAL64) :: CENTX = 0.0D0 , CENTY = 0.0D0, CENTZ = 0.0D0

      CHARACTER(LEN=80) :: INFILE=''
 
! OPEP+HiRE interface
      LOGICAL :: OPEPT = .FALSE., OPEPHIRE_DEBUG = .FALSE.

! Titrations for HiRE
      LOGICAL :: TITRATION = .FALSE.      !Tritration steps?
      INTEGER :: TITMETHOD = 1            !which method
      INTEGER :: TITNSTEPS = 100          !frequency of steps applied
      REAL(KIND = REAl64) :: PH = 7.0     !PH value used for titration
       
! HB for HiRE
      LOGICAL           :: RNAHBT = .FALSE.          !detect h-bonds and save 
      LOGICAL           :: HBSAVET = .FALSE.         !save h-bonds at step 
      INTEGER           :: HBNSTEPS = 1              !frequency for saving h-bonds
! SAXS for HiRE
      LOGICAL           :: SAXST = .FALSE.           !use SAXS
      INTEGER           :: SAXSNSTEPS = 1            !frequency for applying SAXS
      LOGICAL           :: SAXSFORCET = .FALSE.      !compute SASX force at MC transition?
      LOGICAL           :: SAXSMODULT = .FALSE.      !modulate SAXS force with decreasing and periodic factor?
      REAL(KIND=REAL64) :: SAXSINVSIG = 2.5          !broadness of the periodic modulation
      REAL(KIND=REAL64) :: SAXSWAVE = 100.0          !length of the decreasing modulation 
      REAL(KIND=REAL64) :: SAXSOFFI = 0.0D0          !step variable for periodic modulation
      REAL(KIND=REAL64) :: SAXSMODI = 1.0D0          !step variable for decreasing modulation
      LOGICAL           :: SAXSPRINT = .FALSE.       !print curves and scores to unit
      INTEGER           :: SAXSNPRINT = 100          !frequency for printing SAXS
      LOGICAL           :: SAXSSAVET  = .FALSE.      !save computed SAXS at MC transition?
      REAL(KIND=REAL64) :: SAXSMAX = 1.0             !max value for SAXS vector
      LOGICAL           :: SAXSSOLT = .FALSE.        !vacuum (false) or solution (true) calculation
      LOGICAL           :: REFINET = .FALSE.         !refine hydration
      INTEGER           :: NWATLAY = 1               !number of hydration layers
      REAL(KIND=REAL64) :: WATRAD = 2.55             !water radius
      REAL(KIND=REAL64) :: WATW = 0.037              !water contrast


! Dump energies for PiGS
      LOGICAL :: DUMPPIGST = .FALSE.                 !write output for PiGS reading
      REAL(KIND = REAL64) :: EPIGSLIM=HUGE(1.0D0)    !limit for optimisation to save first energy
      REAL(KIND = REAL64) :: EPIGSSAVE(2)            !saved energy
     
! FEBH variables
      LOGICAL           :: FEBHT = .FALSE., SPARSET = .FALSE., SPARSE_BENCH = .FALSE.
      REAL(KIND = REAL64)  :: FEBH_POT_ENE = 0.0D0, FETEMP = 0.0D0, ZERO_THRESH = 0.0D0
      REAL(KIND = REAL64), PARAMETER :: SMALL_DOUBLE = 1.0D-100
      REAL(KIND = REAL64), ALLOCATABLE :: QENERGIES(:), QCOORDINATES(:,:), QPE(:)
      INTEGER           :: FE_FILE_UNIT
      ! Variables for converging to a certain separation of zero and non-zero eigenvalues.
      REAL(KIND = REAL64)  :: MIN_ZERO_SEP = 0.0D0
      INTEGER           :: MAX_ATTEMPTS = 5
      LOGICAL :: USEROT = .FALSE.
      REAL(KIND = REAL64) :: PLANCK
! Needed for partition function calculation
      REAL(KIND = REAL64), ALLOCATABLE :: PFSUM(:)
      REAL(KIND = REAL64), ALLOCATABLE :: EMIN(:), FVIBMIN(:), PFMIN(:), IXMIN(:), IYMIN(:), IZMIN(:), LNFAC(:)

! Rigid body routines
      REAL(KIND = REAL64) :: EPSRIGID=1.0D-3
      LOGICAL :: UPDATERIGIDREFT=.FALSE., HYBRIDMINT=.FALSE., RELAXFQ=.FALSE.
      ! ROTATERIGID variables
      REAL(KIND = REAL64) :: ROTATEFACTOR=1.0D0
      INTEGER :: ROTATERIGIDFREQ=10, ROTRIGIDOFF = 0
      LOGICAL :: ROTATERIGIDT = .FALSE., DOROTATERIGID=.FALSE.
      ! TRANSLATERIGID variables
      REAL(KIND = REAL64) :: TRANSLATEFACTOR = 1D0
      INTEGER :: TRANSLATERIGIDFREQ = 1, TRANSRIGIDOFF = 0
      LOGICAL :: TRANSLATERIGIDT = .FALSE., DOTRANSLATERIGID=.FALSE.

! Freezing routines
      INTEGER :: NFREEZE = 0
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FROZENLIST  !  NATOMS
      LOGICAL :: FREEZE = .FALSE. , FREEZERES = .FALSE., UNFREEZERES = .FALSE.
      LOGICAL :: FREEZEALL = .FALSE.,  UNFREEZEFINALQ  = .FALSE. , FREEZESAVE = .TRUE. 
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FROZENRES  ! NATOMS for safety
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: FROZEN  !  NATOMS

! GROUP ROTATION MOVE PARAMETERS
      INTEGER :: GROUPROTFREQ=1, NGROUPS=0
      LOGICAL :: GROUPROTT=.FALSE., DOGROUPROT=.FALSE., GR_SCALEROT=.FALSE., GR_SCALEPROB=.FALSE.
      LOGICAL :: GROUPROT_SUPPRESS=.FALSE.
      CHARACTER(LEN=10), ALLOCATABLE :: ATOMGROUPNAMES(:)
      CHARACTER(LEN=30) :: GR_SCALEMODE
      INTEGER, ALLOCATABLE :: ATOMGROUPAXIS(:,:)
      REAL(KIND = REAL64), ALLOCATABLE :: ATOMGROUPSCALING(:),ATOMGROUPPSELECT(:)
      LOGICAL, ALLOCATABLE :: ATOMGROUPS(:,:)
      LOGICAL :: SKIPBPT=.FALSE.                         ! skip groups involving paired bases

! Adding pulling force
      LOGICAL :: PULLT = .FALSE.
      INTEGER :: PATOM1 = 0, PATOM2 = 0
      REAL(KIND = REAL64) :: PFORCE = 1.0D1

! Adding harmonic field
      LOGICAL :: HARMONICF = .FALSE.
      REAL(KIND = REAL64) :: HARMONICSTR = 2.5D1
      REAL(KIND = REAL64), ALLOCATABLE :: HARMONICR0(:)
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: HARMONICFLIST  !  NATOMS


! Base-pair inforamtion for large moves
     LOGICAL, ALLOCATABLE, DIMENSION(:,:) :: BP_CURR           ! current base pairs
     TYPE(NUCLEOTIDE), ALLOCATABLE, DIMENSION(:) :: LIST_NUCL  ! list of nucleotides with some additional information
     INTEGER, DIMENSION(2) :: HARMATOMS = (/0,0/)              ! list of atoms for the harmonic field
     INTEGER :: NNUCL                                          ! number of nucleotides
     INTEGER :: NLOOSE = 3                                     ! number of loose nucleotides needed for a loose tail
     INTEGER :: BPHINGEFREQ = 100                              ! frequency of hinge moves
     INTEGER :: NBPHARMOVE = 4                                 ! number of nucleotides from either and that can be chosen for the added spring
     REAL(KIND = REAL64) :: BPTHRESH = 2.3D0                   ! energy threshold for base pairing
     REAL(KIND = REAL64) :: BPDIST = 1.4D0                     ! equilibrium distance for added spring
     REAL(KIND = REAL64) :: BPSTRENGTH = 2.5D1                 ! force constant for added spring
     LOGICAL :: BPHINGET = .FALSE.                             ! logical switch for base pairing moves
     LOGICAL :: NOBPT = .FALSE.                                ! no base pairs logical
     LOGICAL :: HARMONICPOT = .FALSE.                          ! add harmonic field in potential  
     LOGICAL :: DOHINGE = .FALSE.                              ! attempt hinge move in this BH step?
 
! Pull and twist moves
     LOGICAL :: PULLMOVET = .FALSE.                            ! do pulling moves
     LOGICAL :: TWISTMOVET = .FALSE.                           ! do twisting moves
     LOGICAL :: TWISTORPULL = .FALSE.                          ! decide which move was last
     INTEGER :: PULLMFREQ = 50                                 ! frequency of pulling moves
     INTEGER :: TWISTMFREQ = 50                                ! frequency of twisting moves
     INTEGER :: PULLMOFF = 0                                   ! offset subtracted from the number of steps to check freq of moves
     INTEGER :: TWISTMOFF = 0                                  ! offset subtracted from the number of steps to check freq of moves
     LOGICAL :: DOTWIST=.FALSE.                                ! do twisting move in this BH step?
     LOGICAL :: DOPULL=.FALSE.                                 ! do pulling move in this BH step?
     REAL(KIND = REAL64) :: PULLMF = 7.5D0                     ! pulling force applied
     REAL(KIND = REAL64) :: TWISTMF = 7.5D0                    ! twisting force applied

! Harmonic moves
     LOGICAL :: HARMONICMOVET = .FALSE.                        ! do harmonic moves (add spring between bases)
     INTEGER :: HARMOVEFREQ = 50                               ! frequency of moves
     LOGICAL :: DOHARMONIC = .FALSE.                           ! do harmonic move in this BH step?
     REAL(KIND = REAL64) :: HMKF = 2.5D1                       ! spring constant used
     REAL(KIND = REAL64) :: HMDIST = 1.5D0                     ! equilibrium used


! Force adaptation
     LOGICAL :: PADAPTFT = .FALSE.                             !Alter pulling force
     REAL(KIND = REAL64) :: PUPPERF = 2.5D2, PLOWERF = 2.5D0   !Upper and lower limits
     REAL(KIND = REAL64) :: PADAPTSCALE = 1.1D0                !Scaling constant
     LOGICAL :: TADAPTFT = .FALSE.                             !Alter twisting force
     REAL(KIND = REAL64) :: TUPPERF = 2.5D2, TLOWERF = 2.5D0   !Upper and lower limits
     REAL(KIND = REAL64) :: TADAPTSCALE = 1.1D0                !Scaling constant
     
END MODULE COMMONS
