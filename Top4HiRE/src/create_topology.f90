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
         CALL GET_ANGLE_INFO()
         CALL GET_QANGLE_INFO()
         
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
         !angle force constants
         WRITE(TOPUNIT,*) "SECTION ANGLE_FORCE_CONSTANT"
         WRITE(TOPUNIT,'(5F16.8)') AKSPR
         ! angle equilibrium values
         WRITE(TOPUNIT,*) "SECTION ANGLE_EQUIL_VALUE"
         WRITE(TOPUNIT,'(5F16.8)') ATEQ
         ! qangle reference value
         WRITE(TOPUNIT,*) "SECTION QANGLE_REF_ANGLE"
         WRITE(TOPUNIT,'(5F16.8)') QTTS
         ! qangle parameters
         WRITE(TOPUNIT,*) "SECTION QANGLE_A1"
         WRITE(TOPUNIT,'(5F16.8)') QA1
         WRITE(TOPUNIT,*) "SECTION QANGLE_A2"
         WRITE(TOPUNIT,'(5F16.8)') QA2
         WRITE(TOPUNIT,*) "SECTION QANGLE_A3"
         WRITE(TOPUNIT,'(5F16.8)') QA3
         WRITE(TOPUNIT,*) "SECTION QANGLE_A5"
         WRITE(TOPUNIT,'(5F16.8)') QA5
      END SUBROUTINE WRITE_TYPE_DETAILS

      SUBROUTINE WRITE_BONDEDINTS_DETAILS(TOPUNIT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: TOPUNIT
         INTEGER :: I
         !bond information
         WRITE(TOPUNIT,*) "SECTION BONDS"
         WRITE(TOPUNIT,'(12I6)') (BONDS(I,1), BONDS(I,2), BTYPE(I), I=1,NBONDS)
         !angle information
         WRITE(TOPUNIT,*) "SECTION ANGLES"
         WRITE(TOPUNIT,'(12I6)') (ANGLES(I,1), ANGLES(I,2), ANGLES(I,3), ATYPE(I), I=1,NANGLE)     
         !qangle information
         WRITE(TOPUNIT,*) "SECTION QANGLES"
         WRITE(TOPUNIT,'(12I6)') (QANGLES(I,1), QANGLES(I,2), QANGLES(I,3), QTYPE(I), I=1,NQANGLE)               
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
            DO K=1,NBINTRA
               !check whether the type of this bond matches the res type
               IF (.NOT.(BINTRATYPE(K).EQ.RESTYPE(I))) CONTINUE
               !get grain names
               AT1=BINTRA(K)%AT1
               AT2=BINTRA(K)%AT2
               !get their indices
               CALL GET_CG_ID(I,AT1,IDX1)
               CALL GET_CG_ID(I,AT2,IDX2)
               !if either is -1, the bond does not exist
               IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1)) THEN
                  CONTINUE
               ELSE
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
               END IF
            END DO
            !now do interresidue bonds 
            IF (.NOT.TERMINALT) THEN
               DO K=1,NBINTER
                  !check whether the type of this bond matches the res type
                  IF (.NOT.(BINTERTYPE(K).EQ.RESTYPE(I))) CONTINUE
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
                  IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1)) THEN
                     CONTINUE
                  ELSE
                     !otherwise we have a new bond
                     NBONDS = NBONDS + 1
                     DUMMYBONDS(NBONDS,1) = IDX1
                     DUMMYBONDS(NBONDS,2) = IDX2
                     !get the type for the topology
                     IF (TYPEMAP(NBINTRA+K).EQ.-1) THEN
                        !this is a new type
                        NBTYPE = NBTYPE + 1
                        TYPEMAP(NBINTRA+K) = NBTYPE
                     END IF
                     DUMMYTYPE(NBONDS) = TYPEMAP(NBINTRA+K) 
                  END IF
               END DO
            END IF
         END DO
         !now allocate the actual global types and get all bonding information
         ALLOCATE(BONDS(NBONDS,2),BTYPE(NBONDS),BKSPR(NBTYPE),BREQ(NBTYPE))
         BONDS(1:NBONDS,1:2) = 3*(DUMMYBONDS(1:NBONDS,1:2)-1)
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

      SUBROUTINE GET_ANGLE_INFO()
         IMPLICIT NONE
         INTEGER :: I,J,K,IDX1,IDX2, IDX3
         CHARACTER(LEN=4) :: AT1, AT2, AT3
         LOGICAL :: TERMINALT
         INTEGER, ALLOCATABLE :: TYPEMAP(:)
         INTEGER, ALLOCATABLE :: DUMMYTYPE(:), DUMMYANGLE(:,:)
         NANGLE = 0
         ALLOCATE(DUMMYTYPE(NRES*(NAINTRA+NAINTER)),DUMMYANGLE(NRES*(NAINTRA+NAINTER),3),TYPEMAP(NAINTER+NAINTRA))
         !TYPEMAP stores which type each angle is: 1:NAINTRA are the intraresidue angles,
         !NAINTRA+1:NAINTRA+NAINTER are the interresidue types
         !NTYPE stores which type we have right now for new allocations
         !the type mapping is saved for each angle in DUMMYTYPE
         TYPEMAP(1:NAINTER+NAINTRA) = -1
         NATYPE = 0
         DO I=1,NRES
            !check whether this is the end of a chain
            TERMINALT = .FALSE.
            DO J=1,NTERMINI
               IF (TERMINI(J,2).EQ.I) THEN
                  TERMINALT = .TRUE.
                  EXIT
               END IF
            END DO
            !now go over intraresidue angles
            DO K=1,NAINTRA
               !check whether the type of this bond matches the res type
               IF (.NOT.(AINTRATYPE(K).EQ.RESTYPE(I))) CONTINUE
               !get grain names
               AT1=AINTRA(K)%AT1
               AT2=AINTRA(K)%AT2
               AT3=AINTRA(K)%AT3
               !get their indices
               CALL GET_CG_ID(I,AT1,IDX1)
               CALL GET_CG_ID(I,AT2,IDX2)
               CALL GET_CG_ID(I,AT3,IDX3)
               !if either is -1, the bond does not exist
               IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1).OR.(IDX3.EQ.-1)) THEN
                  CONTINUE
               ELSE
                  !otherwise we have a new bond
                  NANGLE = NANGLE + 1
                  DUMMYANGLE(NANGLE,1) = IDX1
                  DUMMYANGLE(NANGLE,2) = IDX2
                  DUMMYANGLE(NANGLE,3) = IDX3
                  !get the type for the topology
                  IF (TYPEMAP(K).EQ.-1) THEN
                     !this is a new type
                     NATYPE = NATYPE + 1
                     TYPEMAP(K) = NATYPE
                  END IF
                  DUMMYTYPE(NANGLE) = TYPEMAP(K) 
               END IF
            END DO
            !now do interresidue bonds 
            IF (.NOT.TERMINALT) THEN
               DO K=1,NAINTER
                  !check whether the type of this bond matches the res type
                  IF (.NOT.(AINTERTYPE(K).EQ.RESTYPE(I))) CONTINUE
                  !get grain names
                  AT1=AINTER(K)%AT1
                  AT2=AINTER(K)%AT2
                  AT3=AINTER(K)%AT3
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
                  IF (INDEX(AT3,"*").GT.0) THEN
                     CALL GET_CG_ID(I+1,AT3,IDX3)
                  ELSE
                     CALL GET_CG_ID(I,AT3,IDX3)
                  END IF
                  !if either is -1, the bond does not exist
                  IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1).OR.(IDX3.EQ.-1)) THEN
                     CONTINUE
                  ELSE
                     !otherwise we have a new bond
                     NANGLE = NANGLE + 1
                     DUMMYANGLE(NANGLE,1) = IDX1
                     DUMMYANGLE(NANGLE,2) = IDX2
                     DUMMYANGLE(NANGLE,3) = IDX3
                     !get the type for the topology
                     IF (TYPEMAP(NAINTRA+K).EQ.-1) THEN
                        !this is a new type
                        NATYPE = NATYPE + 1
                        TYPEMAP(NAINTRA+K) = NATYPE
                     END IF
                     DUMMYTYPE(NANGLE) = TYPEMAP(NAINTRA+K) 
                  END IF
               END DO
            END IF
         END DO
         !now allocate the actual global types and get all bonding information
         ALLOCATE(ANGLES(NANGLE,3),ATYPE(NANGLE),AKSPR(NATYPE),ATEQ(NATYPE))
         ANGLES(1:NANGLE,1:3) = 3*(DUMMYANGLE(1:NANGLE,1:3)-1)
         ATYPE(1:NANGLE) = DUMMYTYPE(1:NANGLE)
         DO I=1,NATYPE
            DO K=1,NAINTRA
               IF (TYPEMAP(K).EQ.I) THEN
                  AKSPR(I) = AINTRA(K)%KSPR
                  ATEQ(I) = AINTRA(K)%TEQ
                  EXIT
               END IF
            END DO
            DO K=1,NAINTER
               IF (TYPEMAP(NAINTRA+K).EQ.I) THEN
                  AKSPR(I) = AINTER(K)%KSPR
                  ATEQ(I) = AINTER(K)%TEQ
                  EXIT  
               END IF             
            END DO
         END DO
      END SUBROUTINE GET_ANGLE_INFO

      SUBROUTINE GET_QANGLE_INFO()
         IMPLICIT NONE
         INTEGER :: I,J,K,IDX1,IDX2, IDX3
         CHARACTER(LEN=4) :: AT1, AT2, AT3
         LOGICAL :: TERMINALT
         INTEGER, ALLOCATABLE :: TYPEMAP(:)
         INTEGER, ALLOCATABLE :: DUMMYTYPE(:), DUMMYQANGLE(:,:)
         NQANGLE = 0
         ALLOCATE(DUMMYTYPE(NRES*(NQINTRA+NQINTER)),DUMMYQANGLE(NRES*(NQINTRA+NQINTER),3),TYPEMAP(NQINTER+NQINTRA))
         !TYPEMAP stores which type each angle is: 1:NQINTRA are the intraresidue angles,
         !NQINTRA+1:NQINTRA+NQINTER are the interresidue types
         !NTYPE stores which type we have right now for new allocations
         !the type mapping is saved for each angle in DUMMYTYPE
         TYPEMAP(1:NQINTER+NQINTRA) = -1
         NQTYPE = 0
         DO I=1,NRES
            !check whether this is the end of a chain
            TERMINALT = .FALSE.
            DO J=1,NTERMINI
               IF (TERMINI(J,2).EQ.I) THEN
                  TERMINALT = .TRUE.
                  EXIT
               END IF
            END DO
            !now go over intraresidue angles
            DO K=1,NQINTRA
               !check whether the type of this bond matches the res type
               IF (.NOT.(QINTRATYPE(K).EQ.RESTYPE(I))) CONTINUE
               !get grain names
               AT1=QINTRA(K)%AT1
               AT2=QINTRA(K)%AT2
               AT3=QINTRA(K)%AT3
               !get their indices
               CALL GET_CG_ID(I,AT1,IDX1)
               CALL GET_CG_ID(I,AT2,IDX2)
               CALL GET_CG_ID(I,AT3,IDX3)
               !if either is -1, the bond does not exist
               IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1).OR.(IDX3.EQ.-1)) THEN
                  CONTINUE
               ELSE
                  !otherwise we have a new bond
                  NQANGLE = NQANGLE + 1
                  DUMMYQANGLE(NQANGLE,1) = IDX1
                  DUMMYQANGLE(NQANGLE,2) = IDX2
                  DUMMYQANGLE(NQANGLE,3) = IDX3
                  !get the type for the topology
                  IF (TYPEMAP(K).EQ.-1) THEN
                     !this is a new type
                     NQTYPE = NQTYPE + 1
                     TYPEMAP(K) = NQTYPE
                  END IF
                  DUMMYTYPE(NQANGLE) = TYPEMAP(K) 
               END IF
            END DO
            !now do interresidue bonds 
            IF (.NOT.TERMINALT) THEN
               DO K=1,NQINTER
                  !check whether the type of this bond matches the res type
                  IF (.NOT.(QINTERTYPE(K).EQ.RESTYPE(I))) CONTINUE
                  !get grain names
                  AT1=QINTER(K)%AT1
                  AT2=QINTER(K)%AT2
                  AT3=QINTER(K)%AT3
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
                  IF (INDEX(AT3,"*").GT.0) THEN
                     CALL GET_CG_ID(I+1,AT3,IDX3)
                  ELSE
                     CALL GET_CG_ID(I,AT3,IDX3)
                  END IF
                  !if either is -1, the bond does not exist
                  IF ((IDX1.EQ.-1).OR.(IDX2.EQ.-1).OR.(IDX3.EQ.-1)) THEN
                     CONTINUE
                  ELSE
                     !otherwise we have a new bond
                     NQANGLE = NQANGLE + 1
                     DUMMYQANGLE(NQANGLE,1) = IDX1
                     DUMMYQANGLE(NQANGLE,2) = IDX2
                     DUMMYQANGLE(NQANGLE,3) = IDX3
                     !get the type for the topology
                     IF (TYPEMAP(NQINTRA+K).EQ.-1) THEN
                        !this is a new type
                        NQTYPE = NQTYPE + 1
                        TYPEMAP(NQINTRA+K) = NQTYPE
                     END IF
                     DUMMYTYPE(NQANGLE) = TYPEMAP(NQINTRA+K) 
                  END IF
               END DO
            END IF
         END DO
         !now allocate the actual global types and get all bonding information
         ALLOCATE(QANGLES(NQANGLE,3),QTYPE(NQANGLE),QTTS(NQTYPE),QA1(NQTYPE), &
                  QA2(NQTYPE),QA3(NQTYPE),QA5(NQTYPE))
         QANGLES(1:NQANGLE,1:3) = 3*(DUMMYQANGLE(1:NQANGLE,1:3)-1)
         QTYPE(1:NQANGLE) = DUMMYTYPE(1:NQANGLE)
         DO I=1,NQTYPE
            DO K=1,NQINTRA
               IF (TYPEMAP(K).EQ.I) THEN
                  QTTS(I) = QINTRA(K)%TTA
                  QA1(I) = QINTRA(K)%A1
                  QA2(I) = QINTRA(K)%A2
                  QA3(I) = QINTRA(K)%A3
                  QA5(I) = QINTRA(K)%A5
                  EXIT
               END IF
            END DO
            DO K=1,NQINTER
               IF (TYPEMAP(NQINTRA+K).EQ.I) THEN
                  QTTS(I) = QINTRA(K)%TTA
                  QA1(I) = QINTRA(K)%A1
                  QA2(I) = QINTRA(K)%A2
                  QA3(I) = QINTRA(K)%A3
                  QA5(I) = QINTRA(K)%A5
                  EXIT  
               END IF             
            END DO
         END DO
      END SUBROUTINE GET_QANGLE_INFO

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
            IF (CGNAMES(I).EQ.SEARCHSTR) THEN
               CGID = I
               EXIT
            END IF
         END DO
      END SUBROUTINE GET_CG_ID

END MODULE CREATE_TOP