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

MODULE PARSE_AA_PDB
   USE TOP_GLOBALS
   USE PREC_HIRE, ONLY: REAL64
   IMPLICIT NONE

   INTEGER :: AANATOMS = 0 ! number of atoms in AA file
   CHARACTER(LEN=4), ALLOCATABLE :: AANAMES(:), AARESNAMES(:)
   REAL(KIND = REAL64), ALLOCATABLE :: XYZAA(:)
   INTEGER, ALLOCATABLE :: AARESSTART(:), AARESFINAL(:)
   INTEGER :: AANTER = 0
   INTEGER, ALLOCATABLE :: AATERMINI(:,:)
   CONTAINS
      SUBROUTINE PARSE_PDB_INPUT(INPUTNAME)
         USE CG_DATA, ONLY: ASSIGN_GRAIN_DATA, GET_RES_DATA
         CHARACTER(LEN=50), INTENT(IN) :: INPUTNAME 
         LOGICAL :: DEBUGT = .FALSE.
         INTEGER :: I
         !get number of atoms, so that we can allocate all necessary arrays
         !within this module we work with all-atom data
         CALL GET_AA_NATOMS(INPUTNAME)
         CALL GET_RESNUM(INPUTNAME)
         !now allocate module variables
         CALL ALLOC_AA_DATA()
         !now parse pdb information properly
         CALL GET_PDB_INFO(INPUTNAME)

         IF (DEBUGT) THEN
            WRITE(*,*) "NATOMS, NRES: ", AANATOMS, NRES
            WRITE(*,*) "AANAMES ", AANAMES
            WRITE(*,*) "NTER, TERMINI: ", AANTER, AATERMINI(1,1), AATERMINI(1,2), AATERMINI(2,1), AATERMINI(2,2) 
            WRITE(*,*) "AARESNAMES: ", AARESNAMES
            WRITE(*,*) "AARESSTART+FINAL ", AARESSTART, AARESFINAL 
         END IF

         !at this point we have all the data, now create the CG data from it
         CALL GET_RES_DATA(AARESNAMES)
         CALL ASSIGN_GRAIN_DATA()
         CALL CREATE_CG_XYZ()
         DO I=1,NATOMS
            WRITE(*,'(A,3F12.7)') CGNAMES(I), XYZCG(3*I-2:3*I)
         END DO
      END SUBROUTINE PARSE_PDB_INPUT

      SUBROUTINE GET_AA_NATOMS(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN
         CHARACTER(LEN=50), INTENT(IN) :: INPUTNAME 
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
         AANATOMS = COUNTER
         WRITE(*,*) " Number of atoms in pdb file: ", AANATOMS
      END SUBROUTINE GET_AA_NATOMS

      SUBROUTINE GET_RESNUM(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN, READLINE
         CHARACTER(LEN=50), INTENT(IN) :: INPUTNAME 
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
 
      SUBROUTINE ALLOC_AA_DATA()
         CALL DEALLOC_AA_DATA()
         ALLOCATE(AANAMES(AANATOMS),XYZAA(3*AANATOMS))
         ALLOCATE(AARESSTART(NRES),AARESFINAL(NRES),AARESNAMES(NRES))
      END SUBROUTINE ALLOC_AA_DATA


      SUBROUTINE DEALLOC_AA_DATA()
         IF (ALLOCATED(AANAMES)) DEALLOCATE(AANAMES)
         IF (ALLOCATED(XYZAA)) DEALLOCATE(XYZAA)
         IF (ALLOCATED(AARESSTART)) DEALLOCATE(AARESSTART)
         IF (ALLOCATED(AARESFINAL)) DEALLOCATE(AARESFINAL)
         IF (ALLOCATED(AARESNAMES)) DEALLOCATE(AARESNAMES)
      END SUBROUTINE DEALLOC_AA_DATA

      SUBROUTINE GET_PDB_INFO(INPUTNAME)
         USE UTILS_IO, ONLY: FILE_OPEN, READLINE
         CHARACTER(LEN=50), INTENT(IN) :: INPUTNAME 
         CHARACTER(LEN=200) :: LINE
         INTEGER :: PDBUNIT, IEND
         LOGICAL :: NEXT
         INTEGER :: CURRATOM, CURRRES, RESID, CURRRESID
         INTEGER :: DUMMYTERMINI(NRES,2)
         INTEGER :: DUMMYNTERMINI

         NEXT=.TRUE.
         CURRATOM = 0
         CURRRES = 0
         CURRRESID = -1
         RESID = 0
         DUMMYNTERMINI = 1
         DUMMYTERMINI(1:NRES,1:2) = -1
         DUMMYTERMINI(DUMMYNTERMINI,1) = 1
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
                  AANAMES(CURRATOM) = LINE(13:16)
                  !xyz coordinates
                  READ(LINE(31:38),'(F8.3)') XYZAA(3*(CURRATOM-1)+1)
                  READ(LINE(39:46),'(F8.3)') XYZAA(3*(CURRATOM-1)+2)
                  READ(LINE(47:54),'(F8.3)') XYZAA(3*(CURRATOM-1)+3)
                  !now deal with the residues
                  !get current id and check whether it is a new residue
                  READ(LINE(23:26),'(I4)') RESID
                  IF (RESID.NE.CURRRESID) THEN
                     CURRRESID = RESID
                     CURRRES = CURRRES + 1
                     IF (CURRRES.GT.1) THEN
                        AARESFINAL(CURRRES-1) = CURRATOM - 1
                     END IF
                     AARESSTART(CURRRES) = CURRATOM
                     AARESNAMES(CURRRES) = LINE(18:20)
                  END IF
               ELSE IF (LINE(1:3).EQ."TER") THEN
                  DUMMYTERMINI(DUMMYNTERMINI,2) = CURRRES
                  IF (CURRRES.LT.NRES) THEN
                     DUMMYNTERMINI = DUMMYNTERMINI + 1
                     DUMMYTERMINI(DUMMYNTERMINI,1) = CURRRES + 1
                  END IF
               END IF
            END IF
         END DO
         CLOSE(PDBUNIT)
         ! add the final res terminus
         AARESFINAL(CURRRES) = CURRATOM
         !allocate terminal arrays
         AANTER = DUMMYNTERMINI
         ALLOCATE(AATERMINI(AANTER,2))
         AATERMINI(1:AANTER,1:2) = DUMMYTERMINI(1:DUMMYNTERMINI,1:2)
      END SUBROUTINE GET_PDB_INFO

      SUBROUTINE CREATE_CG_XYZ()
         IMPLICIT NONE
         INTEGER :: I, J
         INTEGER :: ATOMIDAA, ATOMIDCG

         DO I=1,NRES
            WRITE(*,*) "NEXT RESIDUE: ", I
            IF (RESTYPE(I).EQ.2) THEN !ions only have one set of coordinates
               ATOMIDAA = AARESSTART(J)
               ATOMIDCG = CGSTART(J)
               CALL ASSIGN_AA_TO_CG(ATOMIDAA,ATOMIDCG)
            ELSE IF ((RESTYPE(I).EQ.0).OR.(RESTYPE(I).EQ.1)) THEN !RNA and DNA
               WRITE(*,*) ">>> ", CGRESNAMES(I), CGSTART(I), CGFINAL(I)
               DO J=CGSTART(I),CGFINAL(I) !!!!!!!!!! TODO check definition of CGSTART and CGFINAL
                  WRITE(*,*) CGNAMES(J)
                  IF (CGNAMES(J).EQ."P") THEN
                     CALL GET_ATOMID("P",I,ATOMIDAA)
                     CALL ASSIGN_AA_TO_CG(ATOMIDAA,J)
                  ELSE IF (CGNAMES(J).EQ."O5") THEN
                     CALL GET_ATOMID("O5'",I,ATOMIDAA)
                     CALL ASSIGN_AA_TO_CG(ATOMIDAA,J)
                  ELSE IF (CGNAMES(J).EQ."O3") THEN
                     CALL GET_ATOMID("O3'",I,ATOMIDAA)
                     CALL ASSIGN_AA_TO_CG(ATOMIDAA,J)                     
                  ELSE IF ((CGNAMES(J).EQ."R4").OR.(CGNAMES(J).EQ."S4")) THEN
                     CALL GET_ATOMID("C4'",I,ATOMIDAA)
                     CALL ASSIGN_AA_TO_CG(ATOMIDAA,J)
                  ELSE IF ((CGNAMES(J).EQ."R1").OR.(CGNAMES(J).EQ."S1")) THEN
                     CALL GET_ATOMID("C1'",I,ATOMIDAA)
                     CALL ASSIGN_AA_TO_CG(ATOMIDAA,J)
                  ELSE
                     CALL GET_BASE_COORDS(I,J)
                  END IF
               END DO
            END IF
         END DO
      END SUBROUTINE CREATE_CG_XYZ

      SUBROUTINE GET_BASE_COORDS(RESID,CGID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID,CGID
         CHARACTER(LEN=4), ALLOCATABLE :: GRAINATOMS(:)
         INTEGER :: NATSINGRAIN
         REAL(KIND = REAL64) :: GRAINPOS(3)

         IF (CGNAMES(CGID).EQ."A1") THEN
            NATSINGRAIN = 5
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["C4","C5","N7","C8","N9"]
         ELSE IF (CGNAMES(CGID).EQ."A2") THEN
            NATSINGRAIN = 7
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["N1","C2","N3","C4","C5","C6","N6"]
         ELSE IF (CGNAMES(CGID).EQ."C1") THEN
            NATSINGRAIN = 8
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["N1","C2","O2","N3","N4","C4","C5","C6"]
         ELSE IF (CGNAMES(CGID).EQ."G1") THEN
            NATSINGRAIN = 5
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["C4","C5","N7","C8","N9"]
         ELSE IF (CGNAMES(CGID).EQ."G2") THEN
            NATSINGRAIN = 8
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["N1","C2","N2","N3","C4","C5","C6","O6"]
         ELSE IF (CGNAMES(CGID).EQ."T1") THEN
            NATSINGRAIN = 8
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["N1","C2","N3","C4","C5","C6","O2","O4"]
         ELSE IF (CGNAMES(CGID).EQ."U1") THEN
            NATSINGRAIN = 8
            ALLOCATE(GRAINATOMS(NATSINGRAIN))
            GRAINATOMS = ["N1","C2","N3","C4","C5","C6","O2","O4"]
         END IF
         CALL COMPUTE_GRAINPOS(NATSINGRAIN,GRAINATOMS,RESID,GRAINPOS)
         DEALLOCATE(GRAINATOMS)
         XYZCG(3*(CGID-1)+1:3*(CGID-1)+3) = GRAINPOS(1:3)
      END SUBROUTINE GET_BASE_COORDS

      SUBROUTINE COMPUTE_GRAINPOS(NATSINGRAIN,GRAINATOMS,RESID,GRAINPOS)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATSINGRAIN
         CHARACTER(LEN=4), INTENT(IN) :: GRAINATOMS(NATSINGRAIN)
         INTEGER, INTENT(IN) :: RESID
         REAL(KIND = REAL64), INTENT(OUT) :: GRAINPOS(3)

         INTEGER :: I, ATOMID
         REAL(KIND=REAL64) :: TOTALMASS, ATOMMASS

         GRAINPOS(1:3) = 0.0D0
         TOTALMASS = 0.0D0
         DO I=1,NATSINGRAIN
            CALL GET_ATOMID(GRAINATOMS(I),RESID,ATOMID)
            ATOMMASS = 0.0D0
            IF (INDEX(GRAINATOMS(I),"C").GT.0) THEN
               ATOMMASS = 12.011D0
            ELSE IF (INDEX(GRAINATOMS(I),"N").GT.0) THEN
               ATOMMASS = 14.007D0
            ELSE IF (INDEX(GRAINATOMS(I),"O").GT.0) THEN
               ATOMMASS = 15.999D0
            END IF 
            TOTALMASS = TOTALMASS + ATOMMASS
            GRAINPOS(1) = GRAINPOS(1) + XYZAA(3*(ATOMID-1)+1) * ATOMMASS
            GRAINPOS(2) = GRAINPOS(2) + XYZAA(3*(ATOMID-1)+2) * ATOMMASS
            GRAINPOS(3) = GRAINPOS(3) + XYZAA(3*(ATOMID-1)+3) * ATOMMASS
         END DO
         GRAINPOS(1:3) = GRAINPOS(1:3)/TOTALMASS
      END SUBROUTINE COMPUTE_GRAINPOS

      SUBROUTINE ASSIGN_AA_TO_CG(AAID,CGID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: AAID, CGID
         WRITE(*,*) "AAID, CGID > ",AAID,CGID
         XYZCG(3*(CGID-1)+1) = XYZAA(3*(AAID-1)+1)
         XYZCG(3*(CGID-1)+2) = XYZAA(3*(AAID-1)+2)
         XYZCG(3*(CGID-1)+3) = XYZAA(3*(AAID-1)+3)
      END SUBROUTINE ASSIGN_AA_TO_CG

      !get atom id from name for given residue
      !returns 0 if atom does not exist in residue
      SUBROUTINE GET_ATOMID(ATNAME,RESID,ATOMID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         CHARACTER(*), INTENT(IN) :: ATNAME
         INTEGER, INTENT(OUT) :: ATOMID
         INTEGER :: FIRST, LAST, J1
         WRITE(*,*) " > get atomid> RES: ", RESID, ", atom name: ", ATNAME
         ATOMID = 0
         FIRST=AARESSTART(RESID)
         LAST=AARESFINAL(RESID)
         DO J1=FIRST,LAST
            IF (ADJUSTL(TRIM(AANAMES(J1))).EQ.ADJUSTL(TRIM(ATNAME))) THEN
               ATOMID = J1
               EXIT
            ENDIF
         ENDDO
         RETURN
      END SUBROUTINE GET_ATOMID

END MODULE PARSE_AA_PDB