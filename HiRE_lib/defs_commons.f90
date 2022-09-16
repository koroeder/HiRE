!> @file
!> Contains modules with global variables (PREC_HIRE, VAR_DEFS, VAR_UTILS, SAXS_DEFS, HB_DEFS and NUM_DEFS)

!> Defined precision kinds
MODULE PREC_HIRE
  !> integer 32-bit
  INTEGER, PARAMETER  :: INT32  = SELECTED_INT_KIND(9)
  !> integer 64-bit  
  INTEGER, PARAMETER  :: INT64  = SELECTED_INT_KIND(18)
  !> real, double precision 32-bit
  INTEGER, PARAMETER  :: REAL32 = SELECTED_REAL_KIND(6, 37)
  !> real, double precision 64-bit  
  INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)  
END MODULE PREC_HIRE

!> Global variables for HiRE
MODULE VAR_DEFS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE
  
  !> number of residues
  INTEGER :: NRES
  !> number of CG particles
  INTEGER :: NPARTICLES 
  !> number of degrees of freedom (3*NPARTICLES)
  INTEGER :: NOPT       
  !> number of different particle types
  INTEGER :: NTYPEP     
  !> number of different chains (molecules)
  INTEGER :: NCHAINS    
  
  !> List of charges of all particles
  REAL(KIND = REAL64), ALLOCATABLE :: CHATM(:)  
  !> List of the masses of the particles
  REAL(KIND = REAL64), ALLOCATABLE :: AMASS(:)  
  !> List of particle ids
  INTEGER, ALLOCATABLE :: IAC(:)                
  !> List of indices of first particle in residues
  INTEGER, ALLOCATABLE :: RESSTART(:)           
  !> List of indices of final particle in residues
  INTEGER, ALLOCATABLE :: RESFINAL(:)           
  !> List of types of residue 0-RNA, 1-DNA, 2-protein, 3-other
  INTEGER, ALLOCATABLE :: RESTYPES(:)          
  !> First and last particle for each chain
  INTEGER, ALLOCATABLE :: FRAG_PAR_PTR(:,:)     
  !> First and last residue for each chain 
  INTEGER, ALLOCATABLE :: FRAG_RES_PTR(:,:)     
  !> Number of particles of each residue
  INTEGER, ALLOCATABLE :: RESSIZE(:)            
  !> Number of particles per fragment (molecule)
  INTEGER, ALLOCATABLE :: FRAGSIZE(:)           

  !> List of particle names   
  CHARACTER(4), ALLOCATABLE :: IGRAPH(:)        
  !> List of residue name
  CHARACTER(4), ALLOCATABLE :: RESNAMES(:)      
END MODULE VAR_DEFS

!> Utilities for variables in VAR_DEFS
MODULE VAR_UTILS
  USE PREC_HIRE, ONLY: REAL64
  USE VAR_DEFS 
  IMPLICIT NONE

  CONTAINS 
    !> Allocation for arrays in VAR_DEFS
    SUBROUTINE ALLOC_VARS()
       CALL DEALLOC_VARS()
       ALLOCATE(AMASS(NPARTICLES), CHATM(NPARTICLES), IAC(NPARTICLES), &
                IGRAPH(NPARTICLES), FRAGSIZE(NCHAINS), RESNAMES(NRES), &
                RESSIZE(NRES), RESSTART(NRES), FRAG_PAR_PTR(NCHAINS,2), &
                FRAG_RES_PTR(NCHAINS,2), RESFINAL(NRES), RESTYPES(NRES))
    END SUBROUTINE ALLOC_VARS
    
    !> Deallocation of the arrays in VAR_DEFS
    SUBROUTINE DEALLOC_VARS()
       IF (ALLOCATED(AMASS)) DEALLOCATE(AMASS)
       IF (ALLOCATED(CHATM)) DEALLOCATE(CHATM)
       IF (ALLOCATED(FRAGSIZE)) DEALLOCATE(FRAGSIZE)
       IF (ALLOCATED(FRAG_RES_PTR)) DEALLOCATE(FRAG_RES_PTR)       
       IF (ALLOCATED(FRAG_PAR_PTR)) DEALLOCATE(FRAG_PAR_PTR)       
       IF (ALLOCATED(IAC)) DEALLOCATE(IAC)
       IF (ALLOCATED(IGRAPH)) DEALLOCATE(IGRAPH)       
       IF (ALLOCATED(RESSTART)) DEALLOCATE(RESSTART)
       IF (ALLOCATED(RESFINAL)) DEALLOCATE(RESFINAL)
       IF (ALLOCATED(RESTYPES)) DEALLOCATE(RESTYPES)
       IF (ALLOCATED(RESNAMES)) DEALLOCATE(RESNAMES) 
       IF (ALLOCATED(RESSIZE)) DEALLOCATE(RESSIZE)           
    END SUBROUTINE DEALLOC_VARS    
    
END MODULE VAR_UTILS

!> Module containing global variables for SAXS
MODULE SAXS_DEFS
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE 
   ! saxs serial parameters
   logical :: compute_SAXS_serial, SAXS_save, modulate_SAXS_serial
   real(kind = real64) :: SAXS_invsig
   integer :: saxs_serial_step
   !> output unit
   integer :: SAXSs
   !> output unit
   integer :: SAXSc
   integer :: n_rate_saxs

   !> Boolean deciding if we should compute a SAXS forces next step
   logical :: calc_SAXS_force = .true.
   !> Boolean deciding if we should save SAXS data
   logical :: SAXS_print = .false.
  
   !> use SAXS
   LOGICAL           :: SAXST = .FALSE.           
   !> frequency for applying SAXS
   INTEGER           :: SAXSNSTEPS = 1            
   !> compute SASX force at MC transition?
   LOGICAL           :: SAXSFORCET = .FALSE.      
   !> modulate SAXS force with decreasing and periodic factor?
   LOGICAL           :: SAXSMODULT = .FALSE.      
   !> broadness of the periodic modulation
   REAL(KIND=REAL64) :: SAXSINVSIG = 2.5          
   !> length of the decreasing modulation
   REAL(KIND=REAL64) :: SAXSWAVE = 100.0           
   !> step variable for periodic modulation
   REAL(KIND=REAL64) :: SAXSOFFI = 0.0D0          
   !> step variable for decreasing modulation
   REAL(KIND=REAL64) :: SAXSMODI = 1.0D0          
   !> print curves and scores to unit
   LOGICAL           :: SAXSPRINT = .FALSE.       
   !> frequency for printing SAXS
   INTEGER           :: SAXSNPRINT = 100          
   !> save computed SAXS at MC transition?
   LOGICAL           :: SAXSSAVET  = .FALSE.      
   !> max value for SAXS vector
   REAL(KIND=REAL64) :: SAXSMAX = 1.0             
   !> vacuum (false) or solution (true) calculation
   LOGICAL           :: SAXSSOLT = .FALSE.        
   !> refine hydration
   LOGICAL           :: REFINET = .FALSE.         
   !> number of hydration layers
   INTEGER           :: NWATLAY = 1               
   !> water radius
   REAL(KIND=REAL64) :: WATRAD = 2.55             
   !> water contrast
   REAL(KIND=REAL64) :: WATW = 0.037              
END MODULE SAXS_DEFS

!> Module with the hydrogen bonding global variables
MODULE HB_DEFS
   !> Number of hydrogen-bonds
   INTEGER :: IHB = 0
   !> Is H-bonding initialised?
   LOGICAL :: DO_HB = .FALSE.
   !> Shall the information be saved?
   LOGICAL :: SAVE_HB = .FALSE.
   !> Unit for output file
   INTEGER :: HBDAT
END MODULE HB_DEFS

!> numerical definitions
MODULE NUM_DEFS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE  
  !> value of pi
  REAL(KIND = REAL64), PARAMETER :: PI = 3.141592653589793D0
  !> radian to degree transformation
  REAL(KIND = REAL64), PARAMETER :: RAD2DEG = 57.29577951308232088D0
  !> radian to degree transformation squared
  REAL(KIND = REAL64), PARAMETER :: RAD2 = RAD2DEG * RAD2DEG
END MODULE NUM_DEFS


