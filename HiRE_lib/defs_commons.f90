!Defined precision kinds
MODULE PREC_HIRE
  INTEGER, PARAMETER  :: INT32  = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER  :: INT64  = SELECTED_INT_KIND(18)
  INTEGER, PARAMETER  :: REAL32 = SELECTED_REAL_KIND(6, 37)
  INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)  
END MODULE PREC_HIRE

!Global variables
MODULE VAR_DEFS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE
  
  INTEGER :: NRES       !number of residues
  INTEGER :: NPARTICLES !number of CG particles 
  INTEGER :: NOPT       !number of degrees of freedom (3*NPARTICLES)
  INTEGER :: NTYPEP     !number of different particle types
  INTEGER :: NCHAINS    !number of different particle types
  
  REAL(KIND = REAL64), ALLOCATABLE :: CHATM(:)  !Charge of atoms in simulation
  REAL(KIND = REAL64), ALLOCATABLE :: AMASS(:)  !Masses of the beads
  INTEGER, ALLOCATABLE :: IAC(:)                !Particle ids
  INTEGER, ALLOCATABLE :: RESSTART(:)           !Index of first particle in res
  INTEGER, ALLOCATABLE :: RESFINAL(:)           !Index of final particle in res
  INTEGER, ALLOCATABLE :: RESTYPES(:)           !Type of residue 0-RNA, 1-DNA, 2-protein, 3-other
  INTEGER, ALLOCATABLE :: FRAG_PAR_PTR(:,:)     !First and last particle for each chain
  INTEGER, ALLOCATABLE :: FRAG_RES_PTR(:,:)     !First and last residue for each chain  
  INTEGER, ALLOCATABLE :: RESSIZE(:)            !Number of particles of each residue
  INTEGER, ALLOCATABLE :: FRAGSIZE(:)           !Number of particles per fragment
    
  CHARACTER(4), ALLOCATABLE :: IGRAPH(:)        !Particle names
  CHARACTER(4), ALLOCATABLE :: RESNAMES(:)      !Residue name
  
END MODULE VAR_DEFS

MODULE VAR_UTILS
  USE PREC_HIRE, ONLY: REAL64
  USE VAR_DEFS 
  IMPLICIT NONE

  CONTAINS 
    !Allocation of some generic arrays
    SUBROUTINE ALLOC_VARS()
       CALL DEALLOC_VARS()
       ALLOCATE(AMASS(NPARTICLES), CHATM(NPARTICLES), IAC(NPARTICLES), &
                IGRAPH(NPARTICLES), FRAGSIZE(NCHAINS), RESNAMES(NRES), &
                RESSIZE(NRES), RESSTART(NRES), FRAG_PAR_PTR(NCHAINS,2), &
                FRAG_RES_PTR(NCHAINS,2), RESFINAL(NRES), RESTYPES(NRES))
    END SUBROUTINE ALLOC_VARS
    
    !Deallocation of those arrays
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

MODULE SAXS_DEFS
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE 
   ! saxs serial parameters
   logical :: compute_SAXS_serial, SAXS_save, modulate_SAXS_serial
   real(kind = real64) :: SAXS_invsig
   integer :: saxs_serial_step, SAXSs, SAXSc
   integer :: n_rate_saxs

   ! Booleans deciding if we should compute a SAXS forces next step and save SAXS data
   logical :: calc_SAXS_force = .true.
   logical :: SAXS_print = .false.
  
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
END MODULE SAXS_DEFS

MODULE HB_DEFS
   INTEGER :: IHB = 0
   LOGICAL :: DO_HB = .FALSE.
   LOGICAL :: SAVE_HB = .FALSE.
   INTEGER :: HBDAT
END MODULE HB_DEFS

MODULE NUM_DEFS
  USE PREC_HIRE, ONLY: REAL64
  IMPLICIT NONE  
  
  REAL(KIND = REAL64), PARAMETER :: PI = 3.141592653589793D0
  REAL(KIND = REAL64), PARAMETER :: RAD2DEG = 57.29577951308232088D0
  REAL(KIND = REAL64), PARAMETER :: RAD2 = RAD2DEG * RAD2DEG
END MODULE NUM_DEFS


