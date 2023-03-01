MODULE MD_COMMONS
    USE NUMKIND
    IMPLICIT NONE
    ! logical for REXMD
    LOGICAL :: REXT = .FALSE.
    ! REX method: T- Temperature, H- Hamiltonian
    CHARACTER(LEN=1) :: REXMODE = "H"
    ! number of replicas
    INTEGER :: NREPLICA = 0
    ! Lower bound for replicas
    REAL(KIND = REAL64) :: LOWR = 0.0
    ! Higher bound for replicas
    REAL(KIND = REAL64) :: HIGHR = 1.0
    ! number of tasks
    INTEGER :: NTASKS = 1
    ! task id for MPI
    INTEGER :: TASKID = 0
    ! interval of REX steps
    INTEGER :: NREXSTEPS = 2500
    ! method used in MD, can be VV - Velocity Verlet or LD - Langevin Dynamics
    CHARACTER(LEN=2) :: MDMETHOD = 'VV'
    ! number of MD steps to be taken
    INTEGER :: MDSTEPS = 0
    ! coordinates
    REAL(KIND = REAL64), ALLOCATABLE :: COORDS(:)
    ! accelaration
    REAL(KIND = REAL64), ALLOCATABLE :: ACC(:)
    ! velocity
    REAL(KIND = REAL64), ALLOCATABLE :: VEL(:)
    ! particle masses
    REAL(KIND = REAL64), ALLOCATABLE :: MASSES(:)    
    ! particle names
    CHARACTER(LEN=4), ALLOCATABLE :: ATNAMES(:)
    ! particle else
    CHARACTER(LEN=1), ALLOCATABLE :: ELEMENTS(:)
    ! gamma
    REAL(KIND = REAL64) :: GAMMA = 1.0D-1 
    ! friction parameter
    REAL(KIND = REAL64) :: GFRIC
    ! Langevin scaling parameter, currently not used
    REAL(KIND = REAL64) :: LANGEVINSCALE = 0.1
    ! time step
    REAL(KIND = REAL64) :: DT = 1.0D-2
    ! half a time step
    REAL(KIND = REAL64) :: HDT
    ! temperature
    REAL(KIND = REAL64) :: TEMP = 300.0  
    ! scaling for Hamiltonian
    REAL(KIND= REAL64) :: LAMBDA = 1.0
    ! kinetic energy
    REAL(KIND = REAL64) :: EKIN
    ! potential energy  
    REAL(KIND = REAL64) :: EPOT
    ! number of atoms
    INTEGER :: NATOMS = 0
    ! number of degrees of freedom
    INTEGER :: NOPT = 0
    ! output unit for general output
    INTEGER :: MYUNIT
    ! output unit for structures
    INTEGER :: XUNIT
    ! output unit for energies
    INTEGER :: EUNIT
    ! output unit for RMSD
    INTEGER :: RUNIT
    ! output frequency for coordinates
    INTEGER :: NDUMPX
    ! output frequency for energy
    INTEGER :: NDUMPE   
    ! output frequency for pdb files
    INTEGER :: NDUMPP
    ! output frequency for RMSD
    INTEGER :: NDUMPR
    ! record rmsd?
    LOGICAL :: RMSDT = .FALSE.
    ! Align structures when measuring RMSD?
    LOGICAL :: ALIGNCONFT = .FALSE.
    ! dump pdb files
    LOGICAL :: DUMPPDBT = .FALSE.
    ! Name of topology file        
    CHARACTER(LEN=25) :: TOPNAME = "parameters.top"
    ! Name of scale.dat file
    CHARACTER(LEN=25) :: SCALEDATNAME = "scale_RNA.dat"
    ! Name of initial coordinate file
    CHARACTER(LEN=25) :: COORDSFILE = "start"
    ! Minimise initial structure
    LOGICAL :: MININITIAL = .FALSE.
    ! Use thermalisation
    LOGICAL :: THERMINIT = .FALSE.
    ! Initial temperature
    REAL(KIND = REAL64) ::  TINIT = 1.0D-6
    ! Final temperature
    REAL(KIND = REAL64) ::  TFINAL = 0.616 
    ! Restart simulation from restart file
    LOGICAL :: RESTARTSIMT = .FALSE.
    ! Step number in restart file
    INTEGER :: RESTARTSTEP = 0
    ! Continue simulation at RESTARTSTEP?
    LOGICAL :: CONTINUESIMT = .FALSE.
    ! Restart input file
    CHARACTER(LEN=25) :: RESTARTINPF = "md_restart.dat"
    ! Frequency to dump restart file
    INTEGER :: NDUMPRST = 1000
    SAVE
 END MODULE MD_COMMONS