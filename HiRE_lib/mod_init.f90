!Module to handle initialisation
MODULE MOD_INIT
  USE PREC_HIRE
  IMPLICIT NONE
  CONTAINS
  
    SUBROUTINE TOPOLOGY_INPUT(TOPNAME)
       USE UTILS_IO, ONLY: GETUNIT
       IMPLICIT NONE
       CHARACTER(LEN=*), INTENT(IN) :: TOPNAME !Name of topology file    
       INTEGER :: TOPUNIT
       
       TOPUNIT = GETUNIT()
       OPEN(UNIT=TOPUNIT, FILE=TOPNAME, STATUS='OLD', ACTION='READ')
       CALL READ_TOPOLOGY(TOPUNIT)
       CLOSE(TOPUNIT)
    END SUBROUTINE TOPOLOGY_INPUT
  
    !Read information from topology - assumes file opened and attached to TOPUNIT
    SUBROUTINE READ_TOPOLOGY(TOPUNIT)
       USE UTILS_IO, ONLY: READLINE
       USE VAR_DEFS, ONLY: NRES, NOPT, NPARTICLES, NTYPEP, IGRAPH, RESNAMES, &
                           RESSTART, RESFINAL, FRAG_PAR_PTR, AMASS, IAC, &
                           CHATM, NCHAINS
       USE VAR_UTILS, ONLY: ALLOC_VARS
       USE MOD_BONDS, ONLY: NBONDS, NUMBND, RK, REQ, IB, JB, ICB, ALLOC_BONDS 
       USE MOD_ANGLES, ONLY: NANGLES, NUMANG, TK, TEQ, &
                              IT, JT, KT, ICT, ALLOC_ANGLES 
       USE MOD_DIHEDRALS, ONLY: NDIHS, NPTRA, IP, JP, KP, LP, ICP, PK, PN, &
                              PHASE, ALLOC_DIHS
       
       INTEGER, INTENT(IN) :: TOPUNIT
       INTEGER :: IEND, J
       INTEGER, PARAMETER :: NWORDS = 15
       CHARACTER(LEN=30), DIMENSION(NWORDS) :: WORDSLINE
       CHARACTER(LEN=250) :: THISLINE
       CHARACTER(LEN=30) :: CATNAME
       LOGICAL :: ENDFILET

       !First line is topology name - ignore it
       READ(TOPUNIT, *)
       !Now we iterate, the format is:
       ! SECTION CATNAME
       ! entries ...
       ENDFILET = .FALSE.
       DO WHILE (.NOT.(ENDFILET))
          READ(TOPUNIT, '(A)', IOSTAT=IEND) THISLINE
          IF (IEND.LT.0) THEN
             ENDFILET = .TRUE.
             CYCLE
          ENDIF
          CALL READLINE(THISLINE,NWORDS,WORDSLINE)
          CATNAME = WORDSLINE(2)
          SELECT CASE (CATNAME)
             CASE("DEFINITIONS")
                READ(TOPUNIT,'(12I6)') NPARTICLES, NTYPEP, NBONDS, NANGLES, &
                                       NDIHS, NRES, NUMBND, NUMANG, NPTRA, &
                                       NCHAINS 
                NOPT = 3 * NPARTICLES
                !allocation of arrays from bond, angles and dihedral modules
                CALL ALLOC_BONDS()
                CALL ALLOC_ANGLES()
                CALL ALLOC_DIHS()
                !allocation for general arrays
                CALL ALLOC_VARS()
             CASE("PARTICLE_NAMES")
                READ(TOPUNIT,'(20A4)') (IGRAPH(J), J=1,NPARTICLES)
             CASE("RESIDUE_LABELS")
                READ(TOPUNIT,'(20A4)') (RESNAMES(J), J=1,NRES)
             CASE("RESIDUE_POINTER")
                READ(TOPUNIT,'(12I6)') (RESSTART(J),RESFINAL(J), J=1,NRES)
             CASE("CHAIN_POINTER")
                READ(TOPUNIT,'(12I6)') &
                        (FRAG_PAR_PTR(J,1),FRAG_PAR_PTR(J,2), J=1,NCHAINS)
             CASE("PARTICLE_MASSES")
                READ(TOPUNIT,'(5E16.8)') (AMASS(J), J=1,NPARTICLES)
             CASE("PARTICLE_TYPE")
                READ(TOPUNIT,'(12I6)') (IAC(J), J=1,NPARTICLES)
             CASE("CHARGES")
                READ(TOPUNIT,'(5E16.8)') (CHATM(J), J=1,NPARTICLES)             
             CASE("BOND_FORCE_CONSTANT")
                READ(TOPUNIT,'(5E16.8)') (RK(J), J=1,NUMBND)
             CASE("BOND_EQUIL_VALUE") 
                READ(TOPUNIT,'(5E16.8)') (REQ(J), J=1,NUMBND)
             CASE("ANGLE_FORCE_CONSTANT")
                READ(TOPUNIT,'(5E16.8)') (TK(J), J=1,NUMANG)
             CASE("ANGLE_EQUIL_VALUE")
                READ(TOPUNIT,'(5E16.8)') (TEQ(J), J=1,NUMANG) 
             CASE("DIHEDRAL_FORCE_CONSTANT")
                READ(TOPUNIT,'(5E16.8)') (PK(J), J=1,NPTRA) 
             CASE("DIHEDRAL_PERIODICITY")
                READ(TOPUNIT,'(5E16.8)') (PN(J), J=1,NPTRA) 
             CASE("DIHEDRAL_PHASE")
                READ(TOPUNIT,'(5E16.8)') (PHASE(J), J=1,NPTRA) 
             CASE("BONDS")
                READ(TOPUNIT,'(12I6)') (IB(J), JB(J), ICB(J), J=1,NBONDS)
             CASE("ANGLES")
                READ(TOPUNIT,'(12I6)') &
                           (IT(J), JT(J), KT(J), ICT(J), J=1,NANGLES)
             CASE("DIHEDRALS") 
                READ(TOPUNIT,'(12I6)') &
                           (IP(J), JP(J), KP(J), LP(J), ICP(J), J=1,NDIHS)
          END SELECT
       ENDDO
    END SUBROUTINE READ_TOPOLOGY

    SUBROUTINE READ_SCALE_DAT(SCALEDATNAME)
       USE UTILS_IO, ONLY: GETUNIT
       USE NAPARAMS, ONLY: SCORESIZE, SCORE_RNA 
       IMPLICIT NONE
       
       CHARACTER(LEN=*), INTENT(IN) :: SCALEDATNAME !Name of scale.dat file
       INTEGER :: SCUNIT, I, IDUMMY
       REAL(KIND = REAL64) :: SCORE
       CHARACTER(LEN=100) :: DUMMY
       
       SCUNIT=GETUNIT()
       OPEN(UNIT=SCUNIT, FILE=SCALEDATNAME, STATUS='OLD', ACTION='READ')
       DO I=1,SCORESIZE
          READ(SCUNIT, '(I4,F16.10,A)') IDUMMY, SCORE, DUMMY
          SCORE_RNA(I) = SCORE
       ENDDO
       CLOSE(SCUNIT)   
    END SUBROUTINE READ_SCALE_DAT

    SUBROUTINE INIT_FROM_MODS()
       USE MOD_ANGLES, ONLY: ASSIGN_THETATYPE
       USE MOD_DIHEDRALS, ONLY: INIT_DIHPAR, ASSIGN_PHITYPE
       USE MOD_EXCLV, ONLY: INIT_EXCLV
       USE MOD_DEBYEHUECKEL, ONLY: INIT_DH
       USE MOD_HBONDS, ONLY: SET_HBVARS
       
       CALL INIT_DIHPAR()
       CALL ASSIGN_PHITYPE()
       CALL ASSIGN_THETATYPE()
       CALL INIT_EXCLV()
       CALL INIT_DH()
       CALL SET_HBVARS()

    END SUBROUTINE INIT_FROM_MODS

    SUBROUTINE CREATE_CONSTR()
       USE UTILS_IO, ONLY: FILE_EXIST, GETUNIT 
       USE MOD_RESTRAINTS, ONLY: ALLOC_DISTRESTR, ALLOC_POSRESTR, NRESTS, &
                                 NPOSRES, RESTI, RESTJ, TREST, RESTK, RESTL, &
                                 DRESTL, PRI, PRK, PRX, PRDX
       IMPLICIT NONE
       INTEGER :: RESUNIT, I
       
       IF (FILE_EXIST("restraints.dat")) THEN
          RESUNIT = GETUNIT()
          OPEN(RESUNIT, FILE="restraints.dat", STATUS='OLD', ACTION='READ')
          READ(RESUNIT, '(2I6)') NRESTS, NPOSRES
          CALL ALLOC_DISTRESTR()
          CALL ALLOC_POSRESTR()
          IF (NRESTS.GT.0) THEN
             READ(RESUNIT,*)  !comment line
             READ(RESUNIT, '(2I6,3F15.7,I6)') (RESTI(I), RESTJ(I), RESTK(I), &
                                  RESTL(I), DRESTL(I), TREST(I), I=1,NRESTS)
          ENDIF
          IF (NPOSRES.GT.0) THEN
             READ(RESUNIT,*) !comment line
             READ(RESUNIT, '(I6,7F15.7)') (PRI(I), PRK(I), PRX(1,I), PRX(2,I), &
                        PRX(3,I), PRDX(1,I), PRDX(2,I), PRDX(3,I), I=1,NPOSRES)
          ENDIF
          CLOSE(RESUNIT)
       ELSE
          NRESTS = 0
          NPOSRES = 0
       ENDIF
    END SUBROUTINE CREATE_CONSTR 

    SUBROUTINE INITIALISE_SAXS()
       USE UTILS_IO, ONLY: GETUNIT
       USE SAXS_DEFS, ONLY: SAXSs, SAXSc
       USE SAXS_SCORING, ONLY: SET_SAXS_SCORING

       SAXSs = GETUNIT()
       OPEN(SAXSs, file='SAXS_score.dat', status='unknown',action='write',position='append')
       SAXSc = GETUNIT()
       OPEN(SAXSc, file='SAXS_curve.dat', status='unknown',action='write',position='append')
       CALL SET_SAXS_SCORING()
    END SUBROUTINE INITIALISE_SAXS

END MODULE MOD_INIT 
