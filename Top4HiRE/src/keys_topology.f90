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

!> Module containing global variables
MODULE TOP_GLOBALS
   USE PREC_HIRE
   !> input data mode: one of CG, PDB, SEQ
   CHARACTER(LEN=4) :: MODE
   !> number of CG grains
   INTEGER :: NATOMS = 0
   !> number of nucleotides
   INTEGER :: NRES = 0
   !> list of grain names and residue names
   CHARACTER(LEN=4), ALLOCATABLE :: CGNAMES(:), CGRESNAMES(:)
   !> number of termini
   INTEGER :: NTERMINI = 0
   !> list of termini
   INTEGER, ALLOCATABLE :: TERMINI(:,:)
   !> list of atom types used to look up bonds etc
   INTEGER, ALLOCATABLE :: CGTYPE(:)
   !> list of residue type (RNA, DNA etc.)
   INTEGER, ALLOCATABLE :: RESTYPE(:)
   !> list of start and finish indices for CG residues
   INTEGER, ALLOCATABLE :: CGSTART(:), CGFINAL(:)
   !> list of CG masses
   REAL(KIND = REAL64), ALLOCATABLE :: CGMASS(:)
   !> list of CG charges
   REAL(KIND = REAL64), ALLOCATABLE :: CGCHARGE(:)
   !> coordinates
   REAL(KIND = REAL64), ALLOCATABLE :: XYZCG(:)
   !> bond information
   INTEGER :: NBONDS = 0
   INTEGER :: NBTYPE = 0
   !> bonded atoms
   INTEGER, ALLOCATABLE :: BONDS(:,:)
   !> type for each bond
   INTEGER, ALLOCATABLE :: BTYPE(:)
   !> bond force constants
   REAL(KIND=REAL64), ALLOCATABLE :: BKSPR(:)
   !> bond equilibrium lengths
   REAL(KIND=REAL64), ALLOCATABLE :: BREQ(:)
   !> angle information
   INTEGER :: NANGLE = 0
   INTEGER :: NATYPE = 0
   !> atoms in angle
   INTEGER, ALLOCATABLE :: ANGLES(:,:)
   !> type for each angle
   INTEGER, ALLOCATABLE :: ATYPE(:)
   !> force constants
   REAL(KIND=REAL64), ALLOCATABLE :: AKSPR(:)
   !> equilibrium angle
   REAL(KIND=REAL64), ALLOCATABLE :: ATEQ(:)
   !> qangle information
   INTEGER :: NQANGLE = 0
   INTEGER :: NQTYPE = 0
   !> atoms in qangle
   INTEGER, ALLOCATABLE :: QANGLES(:,:)
   !> type for each qangle
   INTEGER, ALLOCATABLE :: QTYPE(:)
   !> reference angle
   REAL(KIND=REAL64), ALLOCATABLE :: QTTS(:)
   !> potential parameters
   REAL(KIND=REAL64), ALLOCATABLE :: QA1(:), QA2(:), QA3(:), QA5(:)
END MODULE TOP_GLOBALS