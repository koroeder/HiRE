MODULE CREATE_TOP
   USE PREC_HIRE
   USE CG_DATA
   USE TOP_GLOBALS
   USE FF_GLOBALS
   CHARACTER(LEN=40) :: TOPNAME="parameters.top"
   CONTAINS

      SUBROUTINE WRITE_TOPOLOGY()
         USE UTILS_IO, ONLY: GETUNIT
         IMPLICIT NONE
         INTEGER :: TOPUNIT

         CALL GET_BOND_INFO()

         TOPUNIT = GETUNIT()
         OPEN(UNIT=TOPUNIT, FILE=TOPNAME, STATUS='NEW')

         CALL WRITE_PARTICLE_RES_INFO(TOPUNIT)

         CALL WRITE_TYPE_DETAILS(TOPUNIT)

         CALL WRITE_BONDEDINTS_DETAILS(TOPUNIT)
         CLOSE(TOPUNIT)
      END SUBROUTINE WRITE_TOPOLOGY


      SUBROUTINE WRITE_PARTICLE_RES_INFO(TOPUNIT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: TOPUNIT
         INTEGER :: I
         !particle names
         WRITE(TOPUNIT,*) "SECTION PARTICLE_NAMES"
         WRITE(TOPUNIT,'(20A4)') CGNAMES
         !residue names
         WRITE(TOPUNIT,*) "SECTION RESIDUE_NAMES"
         WRITE(TOPUNIT,'(20A4)') CGRESNAMES
         !residue start and finish indices
         WRITE(TOPUNIT,*) "SECTION RESIDUE_POINTER"   
         WRITE(TOPUNIT,'(12I6)') (CGSTART(I),CGFINAL(I), I=1,NRES)
         !termini
         WRITE(TOPUNIT,*) "SECTION CHAIN_POINTER"
         WRITE(TOPUNIT,'(12I6)') (TERMINI(I,1), TERMINI(I,2), I=1,NTERMINI)
         !particle mass
         WRITE(TOPUNIT,*) "SECTION PARTICLE_MASSES"
         WRITE(TOPUNIT,'(5F16.8)') CGMASS
         !particle type
         WRITE(TOPUNIT,*) "SECTION PARTICLE_TYPE"
         WRITE(TOPUNIT,'(12I6)') CGTYPE
      END SUBROUTINE WRITE_PARTICLE_RES_INFO

      SUBROUTINE WRITE_TYPE_DETAILS(TOPUNIT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: TOPUNIT
         INTEGER :: I
         !bond force constants
         WRITE(TOPUNIT,*) "SECTION BOND_FORCE_CONSTANT"
         WRITE(TOPUNIT,'(5F16.8)') BKSPR
         !bond equilibrium distance
         WRITE(TOPUNIT,*) "SECTION BOND_EQUIL_VALUE"
         WRITE(TOPUNIT,'(5F16.8)') BREQ       
      END SUBROUTINE WRITE_TYPE_DETAILS

      SUBROUTINE WRITE_BONDEDINTS_DETAILS(TOPUNIT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: TOPUNIT
         INTEGER :: I
         !bond information
         WRITE(TOPUNIT,*) "SECTION BONDS"
         WRITE(TOPUNIT,'(12I6)') (BONDS(I,1), BONDS(I,2), BTYPE(I), I=1,NBONDS)
      END SUBROUTINE WRITE_BONDEDINTS_DETAILS

      SUBROUTINE GET_BOND_INFO()
         IMPLICIT NONE
         INTEGER :: I,J,K,IDX1,IDX2
         CHARACTER(LEN=4) :: AT1, AT2
         LOGICAL :: TERMINALT
         INTEGER, ALLOCATABLE :: TYPEMAP(:)
         INTEGER, ALLOCATABLE :: DUMMYTYPE(:), DUMMYBONDS(:,:)
         NBONDS = 0
         ALLOCATE(DUMMYTYPE(NRES*(NBINTRA+NBINTER)),DUMMYBONDS(NRES*(NBINTRA+NBINTER),2),TYPEMAP(NBINTER+NBINTRA))
         !TYPEMAP stores which type each bond is 1:NBINTRA are the intraresidue bonds, NBINTRA+1:NBINTRA+NBINTER are the interresidue types
         !NTYPE stores which type we have right now for new allocations
         !the type mappign is saved for each bond in DUMMYTYPE
         TYPEMAP(1:NBINTER+NBINTRA) = -1
         NBTYPE = 0
         DO I=1,NRES
            !check whether this is the end of a chain
            TERMINALT = .FALSE.
            DO J=1,NTERMINI
               IF (TERMINI(J,2).EQ.I) THEN
                  TERMINALT = .TRUE.
                  EXIT
               END IF
            END DO
            !now go over intraresidue bonds
            DO K=1,NBINTRA???????????????
               !check whether the type of this bond matches the res type
               IF (.NOT.(BINTRATYPE(K).EQ.RESTYPE(I))) CONTINUE
               !get grain names
               AT1=BINTRA(K)%AT1
               AT2=BINTRA(K)%AT2
               !get their indices
               CALL GET_CG_ID(I,AT1,IDX1)
               CALL GET_CG_ID(I,AT2,IDX2)
               !if either is -1, the bond does not exist
               IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1)) CONTINUE
               !otherwise we have a new bond
               NBONDS = NBONDS + 1
               DUMMYBONDS(NBONDS,1) = IDX1
               DUMMYBONDS(NBONDS,2) = IDX2
               !get the type for the topology
               IF (TYPEMAP(K).EQ.-1) THEN
                  !this is a new type
                  NBTYPE = NBTYPE + 1
                  TYPEMAP(K) = NBTYPE
               END IF
               DUMMYTYPE(NBONDS) = TYPEMAP(K) 
            END DO
            !now do interresidue bonds 
            IF (.NOT.TERMINALT) THEN
               DO K=1,NBINTER
                  !check whether the type of this bond matches the res type
                  IF (.NOT.(BINTRATYPE(K).EQ.RESTYPE(I))) CONTINUE
                  !get grain names
                  AT1=BINTER(K)%AT1
                  AT2=BINTER(K)%AT2
                  !get their indices - if there is a star in the grain name, it is in the next residue
                  IF (INDEX(AT1,"*").GT.0) THEN
                     CALL GET_CG_ID(I+1,AT1,IDX1)                  
                  ELSE
                     CALL GET_CG_ID(I,AT1,IDX1)
                  END IF
                  IF (INDEX(AT2,"*").GT.0) THEN
                     CALL GET_CG_ID(I+1,AT2,IDX2)
                  ELSE
                     CALL GET_CG_ID(I,AT2,IDX2)
                  END IF                    
                  !if either is -1, the bond does not exist
                  IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1)) CONTINUE
                  !otherwise we have a new bond
                  NBONDS = NBONDS + 1
                  DUMMYBONDS(NBONDS,1) = IDX1
                  DUMMYBONDS(NBONDS,2) = IDX2
                  !get the type for the topology
                  IF (TYPEMAP(NBINTRA+K).EQ.-1) THEN
                     !this is a new type
                     NBTYPE = NBTYPE + 1
                     TYPEMAP(K) = NBTYPE
                  END IF
                  DUMMYTYPE(NBONDS) = TYPEMAP(K) 
               END DO
            END IF
         END DO
         !now allocate the actual global types and get all bonding information
         ALLOCATE(BONDS(NBONDS,2),BTYPE(NBONDS),BKSPR(NBTYPE),BREQ(NBTYPE))
         BONDS(1:NBONDS,1:2) = DUMMYBONDS(1:NBONDS,1:2)
         BTYPE(1:NBONDS) = DUMMYTYPE(1:NBONDS)
         DO I=1,NBTYPE
            DO K=1,NBINTRA
               IF (TYPEMAP(K).EQ.I) THEN
                  BKSPR(I) = BINTRA(K)%KSPR
                  BREQ(I) = BINTRA(K)%REQ
                  EXIT
               END IF
            END DO
            DO K=1,NBINTER
               IF (TYPEMAP(NBINTRA+K).EQ.I) THEN
                  BKSPR(I) = BINTER(K)%KSPR
                  BREQ(I) = BINTER(K)%REQ
                  EXIT  
               END IF             
            END DO
         END DO
      END SUBROUTINE GET_BOND_INFO

      SUBROUTINE GET_CG_ID(RESID,CGNAME,CGID)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: RESID
         CHARACTER(LEN=4), INTENT(IN) :: CGNAME
         INTEGER, INTENT(OUT) :: CGID
         INTEGER :: STARPOS, I
         CHARACTER(LEN=4) :: SEARCHSTR

         CGID = -1
         !check whether the CGname is starred (this is used to indicate which residue atoms belong to)
         STARPOS = INDEX(CGNAME,"*")
         IF (STARPOS.EQ.0) THEN
            SEARCHSTR = CGNAME
         ELSE
            SEARCHSTR = CGNAME
            SEARCHSTR(STARPOS:STARPOS) = " "
         END IF
         DO I=CGSTART(RESID),CGFINAL(RESID)
            IF (CGNAMES(I).EQ.CGNAME) THEN
               CGID = I
               EXIT
            END IF
         END DO
      END SUBROUTINE GET_CG_ID

END MODULE CREATE_TOP