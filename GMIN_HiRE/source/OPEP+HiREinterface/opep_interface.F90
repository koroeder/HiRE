MODULE HIRE_OPEP_INTERFACE_MOD
  
#ifdef __OPEPHIRE
  USE COMMONS, ONLY: DEBUG, MYUNIT, LIST_NUCL, NNUCL, BP_CURR
  USE DEFINITIONS
  USE INPUTmod
  USE CALCFORCES, ONLY: CALCFORCE
  USE PREC_HIRE
#endif
  IMPLICIT NONE
  INTEGER, PARAMETER  :: R64 = SELECTED_REAL_KIND(15, 307)
!******************************************************************************
! Types and parameters for writing pdb ATOM coordinate records to pdb
! files.
  type pdb_atom_data
    character (len = 6)    :: record_name
    integer                :: atom_number
    character (len = 4)    :: atom_name
    character (len = 1)    :: alt_loc_indicator
    character (len = 3)    :: residue_name
    character (len = 1)    :: chain_id
    integer                :: residue_number
    character (len = 1)    :: insertion_code
    real(kind = r64)       :: x_coord
    real(kind = r64)       :: y_coord
    real(kind = r64)       :: z_coord
    real(kind = r64)       :: occupancy
    real(kind = r64)       :: b_factor
    character (len = 2)    :: element
    character (len = 2)    :: charge
  end type pdb_atom_data

! This defines the format to be used when writing ATOM lines for PDB
! files.    
  integer, parameter                   :: pdb_atom_data_size = 15
  type(pdb_atom_data), parameter       :: null_pdb_atom_data = &
    pdb_atom_data('ATOM  ',0,'    ',' ','   ',' ',0,'',0.d0,0.d0,0.d0,1.d0,&
                 &0.d0,'  ','  ')
  character (len=*), parameter         :: atom_string_format = &
  &'(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2,10x,a2,a2)'
!******************************************************************************

  CHARACTER (LEN=4), ALLOCATABLE , SAVE :: RES_NAMES(:), ATOM_NAMES(:)
  INTEGER, ALLOCATABLE ,SAVE :: RES_NUMS(:)
  REAL(KIND = R64), ALLOCATABLE, SAVE :: HiRE_MASSES(:)
#ifdef __OPEPHIRE
  TYPE (t_conformations), SAVE :: CONF
#endif

  CONTAINS

   !kr366> subroutine to determine the number of atoms, identical for both versions  
    SUBROUTINE OPEP_GET_NATOMS(NATOM_)
#ifdef __OPEPHIRE
      USE PORFUNCS
#endif
      INTEGER, INTENT(OUT) :: NATOM_
#ifdef __OPEPHIRE
      !variables needed for creating everything for the output routine later on
      INTEGER :: PDB_FILE_UNIT, GETUNIT, L
      CHARACTER (LEN=10) :: WORD
      CHARACTER (LEN=100) :: LINE
      LOGICAL :: YESNO
#endif
      NATOM_ = 0
#ifdef __OPEPHIRE
      L = 0
      PDB_FILE_UNIT = GETUNIT()
      WRITE(MYUNIT,*) ' HiRE-OPEP_get_atoms> parsing conf_initiale.pdb'
      INQUIRE(FILE='conf_initiale_RNA.pdb',EXIST=YESNO)
      IF (YESNO) THEN
         CALL FILE_OPEN("conf_initiale_RNA.pdb",PDB_FILE_UNIT,.FALSE.)  !open pdb 
      ELSE
         CALL FILE_OPEN("conf_initiale.pdb",PDB_FILE_UNIT,.FALSE.)  !open pdb  
      ENDIF
30    READ(PDB_FILE_UNIT,*,ERR=40,END=40) LINE
      READ(LINE,'(A4)',ERR=30) WORD
      IF (WORD .EQ. "ATOM") THEN
         L = L + 1 
         GOTO 30
      ELSE IF (WORD .EQ. "END") THEN
         GOTO 40
      ELSE
         GOTO 30
      ENDIF
      !close the pdb file and return to main program GMIN
40    CONTINUE
      CLOSE(PDB_FILE_UNIT)
      NATOM_ = L !return number of atoms
      NATOMS = L !set internal number of atoms for OPEP potential
      WRITE(MYUNIT,'(A,I6)') ' HiRE-OPEP_get_atoms> Number of atoms in pdb: ',L
#endif
      RETURN
    END SUBROUTINE OPEP_GET_NATOMS

    SUBROUTINE OPEP_SAXS_INIT()
#ifdef __OPEPHIRE
      CALL INITIALISE_SAXS()
#endif
      RETURN
    END SUBROUTINE OPEP_SAXS_INIT

    SUBROUTINE HIRE_HB_INIT()
#ifdef __OPEPHIRE
      CALL INITIALISE_HB()
#endif
      RETURN
    END SUBROUTINE HIRE_HB_INIT

   !kr366> initiliase potential
    SUBROUTINE OPEP_INIT(NATOM_,COORDS,RNAT,DNAT,PROT)
#ifdef __OPEPHIRE
      USE PORFUNCS
#endif
      INTEGER, INTENT(IN) :: NATOM_
      LOGICAL, INTENT(IN) :: RNAT,DNAT,PROT
      REAL(KIND = R64), INTENT(OUT) :: COORDS(3*NATOM_)
#ifdef __OPEPHIRE
      !variables needed for creating everything for the output routine later on
      INTEGER :: PDB_FILE_UNIT, GETUNIT, L, I, CURRRES, IO_STAT
      CHARACTER (LEN=10) :: WORD,DUMMY
      LOGICAL :: END_PDB, YESNO
#endif
      !dummy initialisation
      COORDS(:) = 0.0
#ifdef __OPEPHIRE
      PDB_FILE_UNIT = GETUNIT()
      !start initialising
      ALLOCATE(RES_NAMES(NATOMS))
      ALLOCATE(RES_NUMS(NATOMS))
      ALLOCATE(ATOM_NAMES(NATOMS))
      ALLOCATE(HiRE_MASSES(NATOMS))
      HiRE_MASSES(:) = 1.0D0
      IF (DEBUG) WRITE(MYUNIT,*) 'HIRE+OPEP_init> call OPEP+HIRE definitions and initialisation'
      CALL INITIALISE(CONF,RNAT,DNAT,PROT) !initialise force field 
      COORDS=pos !pos is a variable from DEFS initiliased in INITIALISE
      ALLOCATE(LIST_NUCL(NNUCL))
      DO I=1,NNUCL
        LIST_NUCL(I)%FATOM = 0
        LIST_NUCL(I)%LATOM = 0
      ENDDO
      ALLOCATE(BP_CURR(NNUCL,NNUCL))
      BP_CURR(:,:) = .FALSE.
      !everything needed for calculations is set up now 
      !before returning to the GMIN routines, initialise inofrmation we need for
      !the output at the end of the run (need residue and atom names)
      IF (DEBUG) WRITE(MYUNIT,*) 'HiRE+OPEP_init> get all information needed for writing output'
      INQUIRE(FILE='conf_initiale_RNA.pdb',EXIST=YESNO)
      IF (YESNO) THEN
         CALL FILE_OPEN("conf_initiale_RNA.pdb",PDB_FILE_UNIT,.FALSE.)  !open pdb 
      ELSE
         CALL FILE_OPEN("conf_initiale.pdb",PDB_FILE_UNIT,.FALSE.)  !open pdb  
      ENDIF
      CURRRES = 0  
30    CALL INPUT(END_PDB, WORD, PDB_FILE_UNIT, .TRUE.)   
      IF (END_PDB) THEN  !reached the end of the pdb, leave the GOTO loop
         GOTO 40
      ENDIF
      IF (WORD .EQ. "ATOM") THEN
         !pdb needs to be correct format as specified below
         CALL READI(L)                       !second entry: atom number
         CALL READA(ATOM_NAMES(L))           !third entry: atom name
         CALL READA(RES_NAMES(L))            !fourth entry: residue name
         IF (RNAT.OR.DNAT) CALL READA(DUMMY) !for RNA we have pdb files with chain identifier, need a smarter way to deal with this
         READ(DUMMY,*,IOSTAT=IO_STAT) RES_NUMS(L)  !lm759> if chain identifier is missing, read res num
         IF (IO_STAT.NE.0) CALL READI(RES_NUMS(L)) !lm759> otherwise, fifth entry: residue number
         IF (CURRRES.LT.RES_NUMS(L)) THEN
            CURRRES = CURRRES + 1
            LIST_NUCL(CURRRES)%FATOM = L
            IF (CURRRES.GT.1) THEN
               LIST_NUCL(CURRRES-1)%LATOM = L - 1
            ENDIF
         END IF
      END IF
      GOTO 30
      !close the pdb file and return to main program
40    CLOSE(PDB_FILE_UNIT)
      LIST_NUCL(NNUCL)%LATOM = NATOM_
      IF (DEBUG) THEN
         WRITE(MYUNIT, *) 'Nucleotide data'
         DO I=1,NNUCL
            WRITE(MYUNIT, *) 'Nucleotide ', I, ', first atom: ', LIST_NUCL(I)%FATOM, ', last atom: ', LIST_NUCL(I)%LATOM
         ENDDO
         WRITE(MYUNIT,*) 'HiRE+OPEP_init> finished potential initialisation'
      ENDIF
#endif
      RETURN
    END SUBROUTINE OPEP_INIT


  SUBROUTINE OPEP_ENERGY_AND_GRADIENT(NATOM_,COORD,GRAD,EREAL)
    INTEGER, INTENT(IN) :: NATOM_
    REAL(KIND = R64), INTENT(IN) :: COORD(3*NATOM_)
    REAL(KIND = R64), INTENT(OUT) :: GRAD(3*NATOM_), EREAL

    GRAD(:) = 0.0D0
    EREAL = 1.0D10
#ifdef __OPEPHIRE
    !call calcforce where we either go for RNA or protein calculation
    !scale factor is set in calcforce as well, use dummy here
    CALL CALCFORCE(1.0D0,COORD,GRAD,EREAL)
    GRAD=-GRAD
#endif
    RETURN
  END SUBROUTINE OPEP_ENERGY_AND_GRADIENT

  SUBROUTINE OPEP_NUM_HESS(NATOM_,COORDS,DELTA,HESSIAN)
    !input and output variable
    INTEGER, INTENT(IN) :: NATOM_
    REAL(KIND = R64), INTENT(IN) :: COORDS(3*NATOM_), DELTA
    REAL(KIND = R64), INTENT(OUT) :: HESSIAN(3*NATOM_,3*NATOM_)
#ifdef __OPEPHIRE
    !variables used internally
    REAL(KIND = R64), DIMENSION(3*NATOMS) :: COORDS_PLUS,COORDS_PLUS2,COORDS_MINUS,COORDS_MINUS2
    REAL(KIND = R64), DIMENSION(3*NATOMS) :: GRAD_PLUS,GRAD_PLUS2,GRAD_MINUS,GRAD_MINUS2
    INTEGER :: I
    REAL(KIND = R64) :: DUMMY_ENERGY
#endif
    HESSIAN(:,:) = 0.0D0
#ifdef __OPEPHIRE  
    !use central finite difference with 4 terms (Richardson interpolation)
    DO I = 1,3*NATOMS
      COORDS_PLUS(:) = COORDS(:)
      COORDS_PLUS(I) = COORDS(I) + DELTA
      CALL OPEP_ENERGY_AND_GRADIENT(NATOMS,COORDS_PLUS,GRAD_PLUS,DUMMY_ENERGY)
      COORDS_PLUS2(:) = COORDS(:)
      COORDS_PLUS2(I) = COORDS(I) + 2.0D0 * DELTA
      CALL OPEP_ENERGY_AND_GRADIENT(NATOMS,COORDS_PLUS2,GRAD_PLUS2,DUMMY_ENERGY)
      COORDS_MINUS(:) = COORDS(:)
      COORDS_MINUS(I) = COORDS(I) - DELTA
      CALL OPEP_ENERGY_AND_GRADIENT(NATOMS,COORDS_MINUS,GRAD_MINUS,DUMMY_ENERGY)
      COORDS_MINUS2(:) = COORDS(:)
      COORDS_MINUS2(I) = COORDS(I) - 2.0D0 * DELTA
      CALL OPEP_ENERGY_AND_GRADIENT(NATOMS,COORDS_MINUS2,GRAD_MINUS2,DUMMY_ENERGY)
      HESSIAN(I,:) = (GRAD_MINUS2(:) - 8.0D0 * GRAD_MINUS(:) + 8.0D0 *GRAD_PLUS(:) - GRAD_PLUS2(:))/(12.0D0*DELTA)
    END DO
#endif
    RETURN
  END SUBROUTINE OPEP_NUM_HESS

  SUBROUTINE OPEP_WRITE_PDB(NATOM_,COORDS,PDB_NAME)
#ifdef __OPEPHIRE
    USE PORFUNCS
#endif
    INTEGER, INTENT(IN) :: NATOM_
    REAL(KIND = R64), INTENT(IN) :: COORDS(3*NATOM_)
    CHARACTER (LEN=*), INTENT(IN) :: PDB_NAME
#ifdef __OPEPHIRE
    INTEGER :: PDB_UNIT, GETUNIT
    INTEGER :: CURR_ATOM
    TYPE(pdb_atom_data) :: CURRENT_ATOM_DATA

    !get unit for pdb file
    PDB_UNIT = GETUNIT()
    OPEN(UNIT=PDB_UNIT, FILE=TRIM(ADJUSTL(PDB_NAME)),ACTION="WRITE")
    WRITE(PDB_UNIT,'(A)') "TITLE"
    DO CURR_ATOM = 1,NATOMS
      CURRENT_ATOM_DATA = null_pdb_atom_data !reset atom data for new atom
      CURRENT_ATOM_DATA % atom_number = CURR_ATOM !atom number
      CURRENT_ATOM_DATA % atom_name = ATOM_NAMES(CURR_ATOM) !atom name
      CURRENT_ATOM_DATA % residue_number = RES_NUMS(CURR_ATOM) !res number
      CURRENT_ATOM_DATA % residue_name = RES_NAMES(CURR_ATOM) !res name
      CURRENT_ATOM_DATA % x_coord = COORDS(3 * CURR_ATOM -2) !x
      CURRENT_ATOM_DATA % y_coord = COORDS(3 * CURR_ATOM -1) !y
      CURRENT_ATOM_DATA % z_coord = COORDS(3 * CURR_ATOM) !z
      WRITE(PDB_UNIT, FMT = atom_string_format) CURRENT_ATOM_DATA
    ENDDO
    WRITE(PDB_UNIT,'(A)') "END"
    CLOSE(PDB_UNIT)
#endif
    RETURN
  END SUBROUTINE OPEP_WRITE_PDB

  SUBROUTINE OPEP_FINISH()
#ifdef __OPEPHIRE
    CALL END_DEFINITIONS()
    DEALLOCATE(RES_NAMES)
    DEALLOCATE(RES_NUMS)
    DEALLOCATE(ATOM_NAMES)
#endif 
    RETURN
  END SUBROUTINE OPEP_FINISH
END MODULE HIRE_OPEP_INTERFACE_MOD 
