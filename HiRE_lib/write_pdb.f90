!> @file
!> Contains PDB_OUT module to produce pdb files

!> Module to handle pdb writing
MODULE PDB_OUT
    USE PREC_HIRE
    IMPLICIT NONE
    !******************************************************************************
    !> Type pdb ATOM records.
    type pdb_atom_data
        !> Record name
        character (len = 6)    :: record_name
        !> Atom number
        integer                :: atom_number
        !> Atom name
        character (len = 4)    :: atom_name
        !> Alternative location indicator
        character (len = 1)    :: alt_loc_indicator
        !> Residue name
        character (len = 3)    :: residue_name
        !> Chain identifier
        character (len = 1)    :: chain_id
        !> Residue number
        integer                :: residue_number
        !> insertion code
        character (len = 1)    :: insertion_code
        !> X coordinate
        real(kind = REAL64)    :: x_coord
        !> Y coordinate
        real(kind = REAL64)    :: y_coord
        !> Z coordinate
        real(kind = REAL64)    :: z_coord
        !> Occupancy
        real(kind = REAL64)    :: occupancy
        !> B factor
        real(kind = REAL64)    :: b_factor
        !> Element
        character (len = 2)    :: element
        !> Charge
        character (len = 2)    :: charge
    end type pdb_atom_data
    !> Number of entries in pdb_atom_data
    integer, parameter                   :: pdb_atom_data_size = 15
    !> Empty entry
    type(pdb_atom_data), parameter       :: null_pdb_atom_data = &
        pdb_atom_data('ATOM  ',0,'    ',' ','   ',' ',0,'',0.d0,0.d0,0.d0,1.d0,&
                 &0.d0,'  ','  ')
    !> Format string for pdb data file
    character (len=*), parameter         :: atom_string_format = &
        &'(a6,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,f8.3,f8.3,f8.3,f6.2,f6.2,10x,a2,a2)'
!******************************************************************************
  
    CONTAINS

        !> Write pdb from coordintes
        !> @brief
        !>
        !> Opens file unit for writing and then calls writepdb subroutine.
        !>
        !> @param[in] COORDS - input coordinates
        !> @param[in] PDBNAME - name of pdb file for output
        !> @param[in] CHAINIDT - should chain identifiers be written?
        SUBROUTINE PDBFROMX(COORDS,PDBNAME,CHAINIDT)
            USE UTILS_IO, ONLY: GETUNIT
            USE VAR_DEFS, ONLY: NPARTICLES
            REAL(KIND=REAL64), INTENT(IN) :: COORDS(3*NPARTICLES)
            CHARACTER(LEN=*), INTENT(IN) :: PDBNAME
            LOGICAL, INTENT(IN) :: CHAINIDT
            INTEGER :: PDBUNIT

            PDBUNIT = GETUNIT()
            OPEN(UNIT=PDBUNIT,FILE=TRIM(ADJUSTL(PDBNAME)),STATUS='UNKNOWN')
            CALL WRITEPDB(COORDS,PDBUNIT,CHAINIDT)
            CLOSE(PDBUNIT)
        END SUBROUTINE PDBFROMX
      
        !> Writes pdb file using the system data loaded in HiRE
        SUBROUTINE WRITEPDB(COORDS,PDB_UNIT,CHAINIDT)
            USE VAR_DEFS, ONLY: NPARTICLES, NCHAINS, RESFINAL, &
                                FRAG_PAR_PTR, IAC, IGRAPH, RESNAMES
            REAL(KIND=REAL64), INTENT(IN) :: COORDS(3*NPARTICLES)
            INTEGER, INTENT(IN) :: PDB_UNIT
            LOGICAL, INTENT(IN) :: CHAINIDT
            INTEGER :: CURR_ATOM, NDUMMY, I, J
            INTEGER :: RES_NUMS(NPARTICLES)
            CHARACTER(LEN=1) :: FRAGID(NPARTICLES)
            TYPE(pdb_atom_data) :: CURRENT_ATOM_DATA
            CHARACTER(LEN=26), PARAMETER :: U = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
            
            ! obtain chain id
            WRITE(*,*) "Create chain ids"
            NDUMMY = 1
            DO I = 1,NPARTICLES
                IF (I.EQ.1) THEN
                    FRAGID(I) = U(NDUMMY:NDUMMY)
                ELSE IF (FRAG_PAR_PTR(NDUMMY,1).EQ.I) THEN
                    NDUMMY = NDUMMY + 1
                    FRAGID(I) = U(NDUMMY:NDUMMY)
                ELSE
                    FRAGID(I) = U(NDUMMY:NDUMMY)
                ENDIF
            ENDDO
            WRITE(*,*) "Done creating chain ids"

            ! obtain residue index for all particles
            NDUMMY = 1
            DO J= 1,NPARTICLES
                RES_NUMS(J) = NDUMMY
                IF (RESFINAL(NDUMMY).EQ.J) THEN
                    NDUMMY = NDUMMY + 1
                ENDIF 
            ENDDO

            NDUMMY = 1
            DO CURR_ATOM = 1,NPARTICLES
                CURRENT_ATOM_DATA = null_pdb_atom_data !reset atom data for new atom
                CURRENT_ATOM_DATA % atom_number = CURR_ATOM !atom number
                CURRENT_ATOM_DATA % atom_name = IGRAPH(CURR_ATOM) !atom name
                CURRENT_ATOM_DATA % residue_number = RES_NUMS(CURR_ATOM) !res number
                CURRENT_ATOM_DATA % residue_name = RESNAMES(RES_NUMS(CURR_ATOM)) !res name
                CURRENT_ATOM_DATA % x_coord = COORDS(3 * CURR_ATOM -2) !x
                CURRENT_ATOM_DATA % y_coord = COORDS(3 * CURR_ATOM -1) !y
                CURRENT_ATOM_DATA % z_coord = COORDS(3 * CURR_ATOM) !z
                IF (CHAINIDT) THEN
                    CURRENT_ATOM_DATA % chain_id = FRAGID(CURR_ATOM)
                ENDIF
                WRITE(PDB_UNIT, FMT = atom_string_format) CURRENT_ATOM_DATA 
                IF ((FRAG_PAR_PTR(NDUMMY,2).EQ.CURR_ATOM).AND.CHAINIDT) THEN
                    NDUMMY = NDUMMY + 1
                    WRITE(PDB_UNIT,'(A)') "TER"
                ENDIF
            ENDDO
            WRITE(PDB_UNIT,'(A)') "END"
        END SUBROUTINE WRITEPDB
     

END MODULE PDB_OUT


