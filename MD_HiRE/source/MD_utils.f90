MODULE MD_COMMONS
    IMPLICIT NONE
    INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)
    ! number of MD steps to be taken
    INTEGER :: MDSTEPS = 0
    ! current MD step
    INTEGER :: CURRSTEP = 0
    ! accelaration
    REAL(KIND = REAL64), ALLOCATABLE :: ACC(:)
    ! velocity
    REAL(KIND = REAL64), ALLOCATABLE :: VEL(:) 
    ! gamma
    REAL(KIND = REAL64) :: GAMMA = 0.0D0 
    ! time step
    REAL(KIND = REAL64) :: DT = 0.0D0
    ! number of atoms
    INTEGER :: NATOMS = 0
    ! output unit for general output
    INTEGER :: MYUNIT
    ! output unit for structures
    INTEGER :: XUNIT
    ! output unit for energies
    INTEGER :: EUNIT
    ! output unit for pdb file
    INTEGER :: PUNIT
    ! output frequency for coordinates
    INTEGER :: NDUMPX
    ! output frequency for energy
    INTEGER :: NDUMPE   
    ! output frequency for pdb files
    INTEGER :: NDUMPP
    ! dump pdb files
    LOGICAL :: DUMPPDBT = .FALSE.
    ! Name of topology file        
    CHARACTER(LEN=25) :: TOPNAME = "parameters.top"
    ! Name of scale.dat file
    CHARACTER(LEN=25) :: SCALEDATNAME = "scale_RNA.dat"
    SAVE
END MODULE MD_COMMONS

MODULE MD_UTILS
   CONTAINS
      SUBROUTINE ALLOC_COMMONS()
         USE MD_COMMONS, ONLY: NATOMS, ACC, VEL
         IMPLICIT NONE
         CALL DEALLOC_COMMONS()
         ALLOCATE(ACC(3*NATOMS))
         ALLOCATE(VEL(3*NATOMS))
      END SUBROUTINE ALLOC_COMMONS

      SUBROUTINE DEALLOC_COMMONS()
        USE MD_COMMONS, ONLY: ACC, VEL
        IMPLICIT NONE        
         IF (ALLOCATED(ACC)) DEALLOCATE(ACC)
         IF (ALLOCATED(VEL)) DEALLOCATE(VEL)
      END SUBROUTINE DEALLOC_COMMONS    
      
      SUBROUTINE REPORT_PARAMS()
         USE MD_COMMONS
         IMPLICIT NONE
      END SUBROUTINE REPORT_PARAMS

END MODULE MD_UTILS