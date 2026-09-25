! Format reference for PDB file
!COLUMNS DATA TYPE FIELD DEFINITION
!-------------------------------------------------------------------------------------
!  1 - 6                  Record name "ATOM "
!  7 - 11   Integer       Atom serial number.
! 13 - 16   Character     Atom name.
!      17   Character     Alternate location indicator.
! 18 - 20   Character     Residue name resName Residue name.
!      22   Character     chainID Chain identifier.
! 23 - 26   Integer       resSeq Residue sequence number.
!      27   AChar         iCode Code for insertion of residues.
! 31 - 38   Real(8.3)     x Orthogonal coordinates for X in Angstroms.
! 39 - 46   Real(8.3)     y Orthogonal coordinates for Y in Angstroms.
! 47 - 54   Real(8.3)     z Orthogonal coordinates for Z in Angstroms.
! 55 - 60   Real(6.2)     occupancy Occupancy.
! 61 - 66   Real(6.2)     tempFactor Temperature factor.
! 77 - 78   LString(2)    element Element symbol, right-justified.
! 79 - 80   LString(2)    charge Charge on the atom.

MODULE PARSE_PDB
   USE TOP_GLOBALS
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE

   INTEGER :: PDBNATOMS = 0 ! number of atoms in AA file
   CHARACTER(LEN=4), ALLOCATABLE :: PDBNAMES(:), PDBRESNAMES(:)
   REAL(KIND = REAL64), ALLOCATABLE :: XYZPDB(:)
   INTEGER, ALLOCATABLE :: AARESSTART(:), PDBRESFINAL(:)
   INTEGER :: PDBNTER = 0
   INTEGER, ALLOCATABLE :: PDBTERMINI(:,:)
   CONTAINS
      SUBROUTINE PARSE_PDB_FILE(INPUTNAME)
         CHARACTER(LEN=256), INTENT(IN) :: INPUTNAME
         LOGICAL :: DEBUGT = .FALSE.
         !get number of atoms, so that we can allocate all necessary arrays
         !within this module we work with all-atom data
         CALL GET_PDB_NATOMS(INPUTNAME)
         CALL GET_RESNUM(INPUTNAME)
         !now allocate module variables
         CALL ALLOC_PDB_DATA()
         !now parse pdb information properly
         CALL GET_PDB_INFO(INPUTNAME)
         !chain starts without phosphate are 5' terminal residues
         CALL MARK_5PRIME_TERMINI()

         IF (DEBUGT) THEN
            WRITE(*,*) "NATOMS, NRES: ", PDBNATOMS, NRES
            WRITE(*,*) "PDBNAMES ", PDBNAMES
            WRITE(*,*) "NTER, TERMINI: ", PDBNTER, PDBTERMINI(1,1), PDBTERMINI(1,2), PDBTERMINI(2,1), PDBTERMINI(2,2)
            WRITE(*,*) "PDBRESNAMES: ", PDBRESNAMES
            WRITE(*,*) "AARESSTART+FINAL ", AARESSTART, PDBRESFINAL
         END IF
      END SUBROUTINE PARSE_PDB_FILE

      SUBROUTINE GET_PDB_NATOMS(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN
         CHARACTER(LEN=256), INTENT(IN) :: INPUTNAME
         CHARACTER(LEN=200) :: LINE
         INTEGER :: COUNTER
         INTEGER :: PDBUNIT, IEND
         LOGICAL :: NEXT
         COUNTER = 0
         NEXT=.TRUE.
         CALL FILE_OPEN(INPUTNAME,PDBUNIT,.FALSE.) !open file read-only
         DO WHILE (NEXT)
            READ(PDBUNIT,'(A)',IOSTAT=IEND) LINE
            IF (IEND.NE.0) THEN
               NEXT=.FALSE.
               CONTINUE
            ELSE
               IF (LINE(1:4).EQ."ATOM") COUNTER = COUNTER + 1
            END IF
         END DO
         CLOSE(PDBUNIT)
         PDBNATOMS = COUNTER
         WRITE(*,*) " Number of atoms in pdb file: ", PDBNATOMS
      END SUBROUTINE GET_PDB_NATOMS

      SUBROUTINE GET_RESNUM(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN, READLINE
         CHARACTER(LEN=256), INTENT(IN) :: INPUTNAME
         CHARACTER(LEN=200) :: LINE
         INTEGER :: COUNTER, CURRENT
         INTEGER :: PDBUNIT, IEND
         LOGICAL :: NEXT
         CHARACTER(LEN=4) :: RESCHAR
         INTEGER :: RESID

         COUNTER = 0
         CURRENT = 0
         NEXT=.TRUE.
         CALL FILE_OPEN(INPUTNAME,PDBUNIT,.FALSE.) !open file read-only
         DO WHILE (NEXT)
            READ(PDBUNIT,'(A)',IOSTAT=IEND) LINE
            IF (IEND.NE.0) THEN
               NEXT=.FALSE.
               CONTINUE
            ELSE
               IF (LINE(1:4).EQ."ATOM") THEN
                  RESCHAR = LINE(23:26)
                  READ(RESCHAR,'(I4)') RESID
                  IF (RESID.NE.CURRENT) THEN
                     CURRENT = RESID
                     COUNTER = COUNTER + 1
                  END IF
               END IF
            END IF
         END DO
         CLOSE(PDBUNIT)
         NRES = COUNTER
      END SUBROUTINE GET_RESNUM

      SUBROUTINE ALLOC_PDB_DATA()
         CALL DEALLOC_PDB_DATA()
         ALLOCATE(PDBNAMES(PDBNATOMS),XYZPDB(3*PDBNATOMS))
         ALLOCATE(AARESSTART(NRES),PDBRESFINAL(NRES),PDBRESNAMES(NRES))
      END SUBROUTINE ALLOC_PDB_DATA

      SUBROUTINE DEALLOC_PDB_DATA()
         IF (ALLOCATED(PDBNAMES)) DEALLOCATE(PDBNAMES)
         IF (ALLOCATED(XYZPDB)) DEALLOCATE(XYZPDB)
         IF (ALLOCATED(AARESSTART)) DEALLOCATE(AARESSTART)
         IF (ALLOCATED(PDBRESFINAL)) DEALLOCATE(PDBRESFINAL)
         IF (ALLOCATED(PDBRESNAMES)) DEALLOCATE(PDBRESNAMES)
      END SUBROUTINE DEALLOC_PDB_DATA

      SUBROUTINE GET_PDB_INFO(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN, READLINE
         CHARACTER(LEN=256), INTENT(IN) :: INPUTNAME
         CHARACTER(LEN=200) :: LINE
         INTEGER :: PDBUNIT, IEND
         LOGICAL :: NEXT
         INTEGER :: CURRATOM, CURRRES, RESID, CURRRESID
         INTEGER :: DUMMYTERMINI(NRES,2)
         INTEGER :: DUMMYNTERMINI
         LOGICAL :: NEWCHAINT
         CHARACTER(LEN=1) :: CURRCHAINID

         NEXT=.TRUE.
         CURRATOM = 0
         CURRRES = 0
         CURRRESID = -1
         RESID = 0
         DUMMYNTERMINI = 0
         DUMMYTERMINI(1:NRES,1:2) = -1
         !the next residue starts a new chain (true for the first residue, and after TER records)
         NEWCHAINT = .TRUE.
         CURRCHAINID = " "
         CALL FILE_OPEN(INPUTNAME,PDBUNIT,.FALSE.) !open file read-only
         DO WHILE (NEXT)
            READ(PDBUNIT,'(A)',IOSTAT=IEND) LINE
            IF (IEND.NE.0) THEN
               NEXT=.FALSE.
               CONTINUE
            ELSE
               IF (LINE(1:4).EQ."ATOM") THEN
                  CURRATOM = CURRATOM + 1
                  !read all information for atom
                  !name
                  PDBNAMES(CURRATOM) = LINE(13:16)
                  !xyz coordinates
                  READ(LINE(31:38),'(F8.3)') XYZPDB(3*(CURRATOM-1)+1)
                  READ(LINE(39:46),'(F8.3)') XYZPDB(3*(CURRATOM-1)+2)
                  READ(LINE(47:54),'(F8.3)') XYZPDB(3*(CURRATOM-1)+3)
                  !now deal with the residues
                  !get current id and check whether it is a new residue
                  READ(LINE(23:26),'(I4)') RESID
                  IF (RESID.NE.CURRRESID) THEN
                     CURRRESID = RESID
                     CURRRES = CURRRES + 1
                     IF (CURRRES.GT.1) THEN
                        PDBRESFINAL(CURRRES-1) = CURRATOM - 1
                        !a change of chain identifier without TER record also starts a new chain
                        IF ((.NOT.NEWCHAINT).AND.(LINE(22:22).NE.CURRCHAINID)) THEN
                           DUMMYTERMINI(DUMMYNTERMINI,2) = CURRRES - 1
                           NEWCHAINT = .TRUE.
                        END IF
                     END IF
                     IF (NEWCHAINT) THEN
                        DUMMYNTERMINI = DUMMYNTERMINI + 1
                        DUMMYTERMINI(DUMMYNTERMINI,1) = CURRRES
                        NEWCHAINT = .FALSE.
                     END IF
                     CURRCHAINID = LINE(22:22)
                     AARESSTART(CURRRES) = CURRATOM
                     PDBRESNAMES(CURRRES) = LINE(18:20)
                  END IF
               ELSE IF (LINE(1:3).EQ."TER") THEN
                  !close the current chain, the next residue (if any) starts a new one
                  IF ((DUMMYNTERMINI.GT.0).AND.(.NOT.NEWCHAINT)) THEN
                     DUMMYTERMINI(DUMMYNTERMINI,2) = CURRRES
                     NEWCHAINT = .TRUE.
                  END IF
               END IF
            END IF
         END DO
         CLOSE(PDBUNIT)
         ! add the final res terminus
         PDBRESFINAL(CURRRES) = CURRATOM
         ! close the last chain if the file does not end with a TER record
         IF (DUMMYNTERMINI.GT.0) THEN
            IF (DUMMYTERMINI(DUMMYNTERMINI,2).EQ.-1) DUMMYTERMINI(DUMMYNTERMINI,2) = CURRRES
         END IF
         !allocate terminal arrays
         PDBNTER = DUMMYNTERMINI
         ALLOCATE(PDBTERMINI(PDBNTER,2))
         PDBTERMINI(1:PDBNTER,1:2) = DUMMYTERMINI(1:DUMMYNTERMINI,1:2)
      END SUBROUTINE GET_PDB_INFO

      !> A nucleotide starting a chain without a phosphate is a 5' terminal residue. Mark it
      !> with a 5 in its name (e.g. G -> G5), so that no P grain is created for it, as if the
      !> 5' terminal name had been used in the input.
      SUBROUTINE MARK_5PRIME_TERMINI()
         INTEGER :: I, J, RESID
         LOGICAL :: HASP
         CHARACTER(LEN=4) :: RNAME

         DO I=1,PDBNTER
            RESID = PDBTERMINI(I,1)
            RNAME = ADJUSTL(PDBRESNAMES(RESID))
            !already named as 5' terminal residue
            IF (INDEX(RNAME,"5").GT.0) CYCLE
            !ions have no phosphate
            IF ((RNAME.EQ."MG").OR.(RNAME.EQ."NA").OR.(RNAME.EQ."K").OR.(RNAME.EQ."CL")) CYCLE
            HASP = .FALSE.
            DO J=AARESSTART(RESID),PDBRESFINAL(RESID)
               IF (ADJUSTL(PDBNAMES(J)).EQ."P") THEN
                  HASP = .TRUE.
                  EXIT
               END IF
            END DO
            IF (.NOT.HASP) THEN
               PDBRESNAMES(RESID) = TRIM(RNAME)//"5"
               WRITE(*,'(A,I6,4A)') "  Residue ", RESID, " (", TRIM(RNAME), &
                                    ") starts a chain without phosphate - treated as 5' terminal ", &
                                    TRIM(PDBRESNAMES(RESID))
            END IF
         END DO
      END SUBROUTINE MARK_5PRIME_TERMINI

END MODULE PARSE_PDB
