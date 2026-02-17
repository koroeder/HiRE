MODULE FF_GLOBALS
   USE PREC_HIRE
   !custom type to store bond information
   TYPE BOND

   END TYPE BOND
   !custom type to store angle information
   TYPE ANGLE

   END TYPE ANGLE
   !custom type to store qangle information
   TYPE QANGLE

   END TYPE QANGLE
   !custom dihedral information type
   TYPE DIHEDRAL

   END TYPE DIHEDRAL
   ! number of bonds between residues
   INTEGER :: NBINTER = 0 
   ! number of bonds within residues
   INTEGER :: NBINTRA = 0
   ! number of angles spanning multiple residues
   INTEGER :: NAINTER = 0 
   ! number of angles within nucleotide
   INTEGER :: NAINTRA = 0
   ! number of qangles spanning multiple nulceotides
   INTEGER :: NQINTER = 0
   ! number of qangles within nucleotide
   INTEGER :: NQINTRA = 0
   ! number of dihedrals spanning multiple nucleotides
   INTEGER :: NDINTER = 0
   ! number of dihedrals within nucleotide
   INTEGER :: NDINTRA = 0
   ! res type for each of the data types above
   INTEGER, ALLOCATABLE :: BINTERTYPE(:), BINTRATYPE(:), AINTERTYPE(:), AINTRATYPE(:), &
                           QINTERTYPE(:), QINTRATYPE(:), DINTERTYPE(:), DINTRATYPE(:)
   ! number of terms for each dihedral
   INTEGER, ALLOCATABLE :: DINTERTERMS(:), DINTRATERMS(:)
   ! maximum number of terms in dihedrals
   INTEGER :: MAXNTERM = 0
   !storing bond information
   TYPE(BOND), ALLOCATABLE :: BINTER(:), BINTRA(:)
   !storing angle information
   TYPE(ANGLE), ALLOCATABLE :: AINTER(:), AINTRA(:)
   !storing qangle information
   TYPE(QANGLE), ALLOCATABLE :: QINTER(:), QINTRA(:)
   !storing sihedrals
   TYPE(DIHEDRAL), ALLOCATABLE :: DINTER(:,:), DINTRA(:,:)


END MODULE FF_GLOBALS