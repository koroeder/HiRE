MODULE MD_COMMONS
   USE NUMKIND
   IMPLICIT NONE
   ! number of MD steps to be taken
   INTEGER :: MDSTEPS = 0
   ! coordinates
   REAL(KIND = REAL64), ALLOCATABLE :: COORDS(:)
   ! accelaration
   REAL(KIND = REAL64), ALLOCATABLE :: ACC(:)
   ! velocity
   REAL(KIND = REAL64), ALLOCATABLE :: VEL(:)
   ! atomic masses
   REAL(KIND = REAL64), ALLOCATABLE :: MASSES(:)    
   ! gamma
   REAL(KIND = REAL64) :: GAMMA = 1.0D-1 
   ! friction parameter
   REAL(KIND = REAL64) :: GFRIC
   ! time step
   REAL(KIND = REAL64) :: DT = 1.0D-2
   ! half a time step
   REAL(KIND = REAL64) :: HDT
   ! temperature
   REAL(KIND = REAL64) :: TEMP = 300.0  
   ! kinetic energy
   REAL(KIND = REAL64) :: EKIN
   ! potential energy  
   REAL(KIND = REAL64) :: EPOT
   ! number of atoms
   INTEGER :: NATOMS = 0
   ! output unit for general output
   INTEGER :: MYUNIT
   ! output unit for structures
   INTEGER :: XUNIT
   ! output unit for energies
   INTEGER :: EUNIT
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
   ! Name of initial coordinate file
   CHARACTER(LEN=25) :: COORDSFILE = "start"
   SAVE
END MODULE MD_COMMONS

MODULE MD_UTILS
   CONTAINS
      SUBROUTINE ALLOC_COMMONS()
         USE MD_COMMONS, ONLY: NATOMS, ACC, VEL, MASSES, COORDS
         IMPLICIT NONE
         CALL DEALLOC_COMMONS()
         ALLOCATE(COORDS(3*NATOMS))
         ALLOCATE(ACC(3*NATOMS))
         ALLOCATE(VEL(3*NATOMS))
         ALLOCATE(MASSES(NATOMS))
      END SUBROUTINE ALLOC_COMMONS

      SUBROUTINE DEALLOC_COMMONS()
         USE MD_COMMONS, ONLY: ACC, VEL, COORDS, MASSES
         IMPLICIT NONE        
         IF (ALLOCATED(COORDS)) DEALLOCATE(COORDS)
         IF (ALLOCATED(ACC)) DEALLOCATE(ACC)
         IF (ALLOCATED(VEL)) DEALLOCATE(VEL)
         IF (ALLOCATED(MASSES)) DEALLOCATE(MASSES)
      END SUBROUTINE DEALLOC_COMMONS    
      
      SUBROUTINE REPORT_PARAMS()
         USE MD_COMMONS
         IMPLICIT NONE
         WRITE(MYUNIT,'(A)') " Molecular dynamics simulation for HiRE "
         WRITE(MYUNIT,'(A)') " ______________________________________ "
         WRITE(MYUNIT,'(A)') " "
         WRITE(MYUNIT,'(A,I10,A)') " settings> Run MD simulation for ", MDSTEPS, "steps"
         WRITE(MYUNIT,'(A,F6.2,A)') " settings> Time step for simulation: ", DT, "fs"
         WRITE(MYUNIT,'(A,F6.2)') " settings> Gamma value for Langevin dynamics: ", GAMMA
         WRITE(MYUNIT,'(A,F8.2)') " settings> Temperature for MD simulation: ", TEMP
         WRITE(MYUNIT,'(A)') " "
      END SUBROUTINE REPORT_PARAMS

      SUBROUTINE MD_START()
         USE MD_COMMONS, ONLY: MYUNIT
         USE FILE_UTILS, ONLY: FILE_OPEN, FILE_EXIST

         IF (FILE_EXIST("mddata")) THEN
            CALL FILE_OPEN("mdout.log",MYUNIT,.TRUE.)
         ELSE
            WRITE(*,'(A)') " Cannot locate input file"
            STOP
         END IF
      END SUBROUTINE MD_START

      SUBROUTINE SET_DERIVED_PARAMS()
         USE MD_COMMONS, ONLY: DT, HDT, GAMMA, GFRIC
         IMPLICIT NONE
         HDT = 0.5*DT
         GFRIC = 1.0D0 - GAMMA*HDT
      END SUBROUTINE SET_DERIVED_PARAMS
END MODULE MD_UTILS