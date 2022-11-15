MODULE BP_MOVES_MOD
  USE COMMONS, ONLY: NATOMS, DEBUG, MYUNIT, BP_CURR, LIST_NUCL, NOBPT, NNUCL, NLOOSE
  USE PREC
  USE VEC3
  USE ROTATIONS, ONLY: ROT_AA2MX
  USE HIRE_INTERFACE, ONLY: DUMP_PDB

  IMPLICIT NONE
  INTEGER, SAVE :: BP_STEPS = 0

  CONTAINS

  SUBROUTINE ANALYSE_BP(LFIVET,LTHREET,NUCID5,NUCID3)
    INTEGER, INTENT(OUT) :: NUCID5, NUCID3        !ID of furthest loose nucleotide from 5' and 3'
    LOGICAL, INTENT(OUT) :: LFIVET,LTHREET        !Loose 5' or 3' end (at least NLOOSE nucleotides)
    INTEGER :: I1,I2

    LFIVET = .FALSE.
    LTHREET = .FALSE.
    NUCID5 = 0
    NUCID3 = 0
    !find length of loose end starting from 5' end
    DO I1=1,NNUCL-1      
       DO I2=I1+1,NNUCL
          IF (BP_CURR(I1,I2)) THEN
             NUCID5 = I1-1
             WRITE(MYUNIT,'(A,I8,A,I8)') ' bpmoves> First base pair from 5-end between ', I1, ' and ', I2
             GOTO 20
          ENDIF
       ENDDO
    ENDDO
    !at this stage we only leave the loop normally if NUCID5 is 0, i.e. there is no base pair
    !in this case, we exit and set the no base pair variable accordingly
    IF (NUCID5.EQ.0) THEN
        NOBPT = .TRUE.
        RETURN
    END IF
20  IF (NUCID5.GE.NLOOSE) LFIVET=.TRUE.
    !find length of loose end starting from 3' end
    DO I1=NNUCL,2,-1    
       DO I2=NNUCL-1,1,-1
          IF (BP_CURR(I1,I2)) THEN
             NUCID3 = I1+1
             WRITE(MYUNIT,'(A,I8,A,I8)') ' bpmoves> First base pair from 3-end between ', I1, ' and ', I2
             GOTO 30
          ENDIF
       ENDDO
    ENDDO
30  IF ((NNUCL-NUCID3).GE.NLOOSE) LTHREET=.TRUE.
    RETURN
  END SUBROUTINE ANALYSE_BP
  
  SUBROUTINE  MOVE_TAILS(X,ITER)
    USE COMMONS
    REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
    INTEGER, INTENT(IN) :: ITER
    REAL(KIND=REAL64) :: CDUMMY(3*NATOMS), EDUMMY, RANDOM, DPRAND
    INTEGER :: U, V
    LOGICAL :: CFLAG

    !call random number generator once to initiate if it wasn't seeded correctly
    RANDOM = DPRAND()

    ! select a random number from 0 to 1:
    RANDOM = DPRAND()
    ! now generate a random integer U between 1 and NLOOSE
    U = FLOOR(NBPHARMOVE*RANDOM+1)
    ! repeat the same for V, used for the 3 tail
    RANDOM = DPRAND()
    V = NNUCL - FLOOR(NBPHARMOVE*RANDOM)
    ! we add a harmonic potential between the atoms, and minimise
    HARMATOMS(1) = LIST_NUCL(U)%LATOM
    HARMATOMS(2) = LIST_NUCL(V)%LATOM
    WRITE(MYUNIT, '(A,I8,A,I8)') " bpmoves> Add spring between nucleotides ", U , " and " , V
    HARMONICPOT = .TRUE.
    CDUMMY(:) = X(:)
    CALL MYLBFGS(3*NATOMS,MUPDATE,CDUMMY,.FALSE.,BQMAX,CFLAG,EDUMMY,MAXIT,ITER,.TRUE.)
    IF (CFLAG) THEN !if converged do some tests
       WRITE(MYUNIT, '(A)') " bpmoves> Converged minimisation, remove spring"
       X(:) = CDUMMY(:)
    ELSE
       ! attept minimisation again with more steps
       CALL MYLBFGS(3*NATOMS,MUPDATE,CDUMMY,.FALSE.,BQMAX,CFLAG,EDUMMY,2*MAXIT,ITER,.FALSE.)
       IF (CFLAG) THEN
          WRITE(MYUNIT, '(A)') " bpmoves> Converged minimisation, remove spring"
          X(:) = CDUMMY(:)
       ELSE
          WRITE(MYUNIT, '(A)') " bpmoves> Could not converge structure, no move attempted"
       ENDIF
    ENDIF
    HARMONICPOT = .FALSE.
    RETURN
  END SUBROUTINE

  SUBROUTINE MOVE_BOTH_TAILS(X, NUCIDF, NUCIDT)
    INTEGER, INTENT(IN) :: NUCIDF, NUCIDT
    REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
    INTEGER :: NATS, NTAIL
    REAL(KIND = REAL64), ALLOCATABLE :: XTEMP(:)

    ! start with five tail
    NATS = LIST_NUCL(NUCIDT-1)%LATOM ! atoms in 3' tail removed
    ALLOCATE(XTEMP(3*NATS))
    XTEMP(:) = X(1:3*NATS)
    CALL MOVE_FTAIL(XTEMP, NATS, NUCIDF, .FALSE.)
    X(1:3*NATS) = XTEMP(:)
    DEALLOCATE(XTEMP)
    ! then continue with the three tail
    NTAIL = LIST_NUCL(NUCIDF)%LATOM
    NATS = NATOMS - NTAIL
    ALLOCATE(XTEMP(3*NATS))
    XTEMP(:) = X(3*NTAIL+1:3*NATOMS)
    CALL MOVE_TTAIL(XTEMP, NATS, NUCIDT, .FALSE.)
    X(3*NTAIL+1:3*NATOMS) = XTEMP(:)
    DEALLOCATE(XTEMP)

  END SUBROUTINE MOVE_BOTH_TAILS

  SUBROUTINE MOVE_FTAIL(X, NATS, NUCID, PRINTPDB)
    INTEGER, INTENT(IN) :: NUCID, NATS
    REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATS)
    LOGICAL, INTENT(IN) :: PRINTPDB
    REAL(KIND = REAL64), ALLOCATABLE :: XTEMP(:)
    REAL(KIND = REAL64) :: BPVEC(3), BBVEC1(3), BBVEC2(3), BBPERP1(3), BBPERP2(3)
    REAL(KIND = REAL64) :: ROTVEC(3), ROTMAT(3,3), ANGLE, CENTRE(3), TARGET2(3), BBTAR(3)
    REAL(KIND = REAL64) :: COM(3)
    INTEGER :: I, J, NUC1AT, NUC2AT, NATTEMP, FIRSTATOM, BPID, CENTREID
    REAL(KIND = REAL64), PARAMETER :: ZERO = 0.0D0, SCALEDIST = 1.2D0
    CHARACTER(LEN=6)    :: J1_STRING

  !1. Bring orientation between first and second loose nucleotide approx. parallel to last base pair
    !First find the vector parallel to the last base pair (approximately)
    BPID = 0
    BPVEC(:) = ZERO
    DO I=NUCID+1,NNUCL
       IF (BP_CURR(NUCID+1,I)) THEN
          NUC1AT = LIST_NUCL(I)%LATOM     
          NUC2AT = LIST_NUCL(NUCID+1)%LATOM
          BPID=I
          BPVEC(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
          BPVEC(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
          BPVEC(3) = X(3*NUC2AT) - X(3*NUC1AT)
          EXIT
       ENDIF
    ENDDO
    IF (VEC_LEN(BPVEC).EQ.0.0D0) THEN
       WRITE(MYUNIT, '(A)') " bpmoves> vector for last base pair is zero, no move attempted"
       RETURN
    END IF
    BP_STEPS = BP_STEPS + 1
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_a."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

    !Now find vector between nucid and nucid - 1 (BBVEC1), and vector between nucid +1 and nucid (BBVEC2)
    BBVEC1(:) = ZERO
    BBVEC2(:) = ZERO
    BBPERP1(:) = ZERO
    BBPERP2(:) = ZERO
    NUC1AT = LIST_NUCL(NUCID)%FATOM     
    NUC2AT = LIST_NUCL(NUCID-1)%FATOM    
    BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
    NUC1AT = LIST_NUCL(NUCID+1)%FATOM     
    NUC2AT = LIST_NUCL(NUCID)%FATOM    
    BBVEC2(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC2(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC2(3) = X(3*NUC2AT) - X(3*NUC1AT)
    !Perpendicular componennt for BBVEC1
    CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
    !Perpendicular component for BPVEC
    CALL PERPVEC(BPVEC,BBVEC2,BBPERP2)
    ! Angle between these two perpendicular components is the rotational angle for the first rotation
    ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
    ROTVEC = ANGLE*VEC_NORM(BBVEC2)
    ! Create rotational matrix
    ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
    ! Select subset of coordinates
    NATTEMP = LIST_NUCL(NUCID)%LATOM !number of atoms in XTEMP
    ALLOCATE(XTEMP(3*NATTEMP))
    XTEMP(:) = X(1:3*NATTEMP)
    CENTREID = LIST_NUCL(NUCID)%FATOM
    ! find new centre
    CENTRE(:) = XTEMP(3*CENTREID-2:3*CENTREID)
    ! Apply rotation
    DO J=1,NATTEMP
       XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
       XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
       XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
    ENDDO
    ! Apply update of coords to X
    X(1:3*NATTEMP) = XTEMP(:)
    DEALLOCATE(XTEMP)
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_b."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

  !2. Rotate third free nucleotide to be approximately in line with base pair
    ! define target coordinate for third loose nucleotide
    NUC1AT = LIST_NUCL(NUCID+1)%FATOM
    TARGET2(1) = X(3*NUC1AT-2) + SCALEDIST*BPVEC(1)
    TARGET2(2) = X(3*NUC1AT-1) + SCALEDIST*BPVEC(2)
    TARGET2(3) = X(3*NUC1AT) + SCALEDIST*BPVEC(3)      
    !Now find vector between nucid - 1 and nucid - 2 (BBVEC1), and vector between nucid and nucid -1 (BBVEC2)
    BBVEC1(:) = ZERO
    BBVEC2(:) = ZERO
    BBPERP1(:) = ZERO
    BBPERP2(:) = ZERO
    NUC1AT = LIST_NUCL(NUCID-1)%FATOM     
    NUC2AT = LIST_NUCL(NUCID-2)%FATOM    
    BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
    NUC1AT = LIST_NUCL(NUCID)%FATOM     
    NUC2AT = LIST_NUCL(NUCID-1)%FATOM    
    BBVEC2(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC2(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC2(3) = X(3*NUC2AT) - X(3*NUC1AT)
    ! set vector to TARGET
    NUC1AT = LIST_NUCL(NUCID-1)%FATOM 
    BBTAR(1) = TARGET2(1) - X(3*NUC1AT-2) 
    BBTAR(2) = TARGET2(2) - X(3*NUC1AT-1)
    BBTAR(3) = TARGET2(3) - X(3*NUC1AT)
    !Perpendicular componennt for BBVEC1
    CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
    !Perpendicular component for BBTAR
    CALL PERPVEC(BBTAR,BBVEC2,BBPERP2)
    ! Angle between these two perpendicular components is the rotational angle for the second rotation
    ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
    ROTVEC = ANGLE*VEC_NORM(BBVEC2)
    ! Create rotational matrix
    ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
    ! Select subset of coordinates
    NATTEMP = LIST_NUCL(NUCID-1)%LATOM !number of atoms in XTEMP
    ALLOCATE(XTEMP(3*NATTEMP))
    XTEMP(:) = X(1:3*NATTEMP)
    CENTREID = LIST_NUCL(NUCID-1)%FATOM
    ! find new centre
    CENTRE(:) = XTEMP(3*CENTREID-2:3*CENTREID)
    ! Apply rotation
    DO J=1,NATTEMP
       XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
       XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
       XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
    ENDDO
    ! Apply update of coords to X
    X(1:3*NATTEMP) = XTEMP(:)
    DEALLOCATE(XTEMP)
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_c."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

  !3. If we have more than 3 loose nucleotides, rotate the next bond closer to the CoM of rest of the molecule
    IF ((NUCID-3).GT.0) THEN
       ! first determine CoM of the rest of the molecule
       FIRSTATOM = LIST_NUCL(NUCID+1)%FATOM
       ALLOCATE(XTEMP(3*(NATS-FIRSTATOM+1)))
       XTEMP(:) = X(3*FIRSTATOM-2:3*NATS)
       COM(:) = ZERO
       DO I=1,NATS-FIRSTATOM+1
          COM(1) = COM(1) + XTEMP(3*I-2)
          COM(2) = COM(2) + XTEMP(3*I-1)
          COM(3) = COM(3) + XTEMP(3*I)
       ENDDO
       COM(:) = COM(:)/(NATS-FIRSTATOM+1)
       DEALLOCATE(XTEMP)
       ! target: vector parallel to middle of last base pair to COM
       NUC1AT = LIST_NUCL(BPID)%LATOM     
       NUC2AT = LIST_NUCL(NUCID+1)%LATOM  
       BPVEC(1) = COM(1) - (X(3*NUC2AT-2) + X(3*NUC1AT-2))/2
       BPVEC(2) = COM(2) - (X(3*NUC2AT-1) + X(3*NUC1AT-1))/2
       BPVEC(3) = COM(3) - (X(3*NUC2AT) + X(3*NUC1AT))/2
       !Now find vector between nucid - 2 and nucid - 3 (BBVEC1), and vector BBVEC1 and BPVEC (BBVEC2)
       BBVEC1(:) = ZERO
       BBVEC2(:) = ZERO
       BBPERP1(:) = ZERO
       BBPERP2(:) = ZERO
       NUC1AT = LIST_NUCL(NUCID-2)%FATOM     
       NUC2AT = LIST_NUCL(NUCID-3)%FATOM    
       BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
       BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
       BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
       !rotational reference is between BBVEC1 and TARGET
       BBVEC2(1) = BBVEC1(1) + BPVEC(1)
       BBVEC2(2) = BBVEC1(2) + BPVEC(2)
       BBVEC2(3) = BBVEC1(3) + BPVEC(3)
       !Perpendicular componennt for BBVEC1
       CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
       !Perpendicular component for BBTAR
       CALL PERPVEC(BPVEC,BBVEC2,BBPERP2)
       ! Angle between these two perpendicular components is the rotational angle for the first rotation
       ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
       ROTVEC = ANGLE*VEC_NORM(BBVEC2)
       ! Create rotational matrix
       ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
       ! Select subset of coordinates
       NATTEMP = LIST_NUCL(NUCID-2)%LATOM !number of atoms in XTEMP
       ALLOCATE(XTEMP(3*NATTEMP))
       XTEMP(:) = X(1:3*NATTEMP)
       ! find new centre
       CENTRE(:) = XTEMP(1:3)
       ! Apply rotation
       DO J=1,NATTEMP
          XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
          XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
          XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
       ENDDO
       ! Apply update of coords to X
       X(1:3*NATTEMP) = XTEMP(:)
       DEALLOCATE(XTEMP)
    ENDIF
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_d."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF
    RETURN
  END SUBROUTINE MOVE_FTAIL

  SUBROUTINE MOVE_TTAIL(X, NATS, NUCID, PRINTPDB)
    INTEGER, INTENT(IN) :: NUCID, NATS
    REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATS)
    LOGICAL, INTENT(IN) :: PRINTPDB
    REAL(KIND = REAL64), ALLOCATABLE :: XTEMP(:)
    REAL(KIND = REAL64) :: BPVEC(3), BBVEC1(3), BBVEC2(3), BBPERP1(3), BBPERP2(3)
    REAL(KIND = REAL64) :: ROTVEC(3), ROTMAT(3,3), ANGLE, CENTRE(3), TARGET2(3), BBTAR(3)
    REAL(KIND = REAL64) :: COM(3)
    INTEGER :: I, J, NUC1AT, NUC2AT, NATTEMP, BPID, FIRSTAT, LASTATOM
    REAL(KIND = REAL64), PARAMETER :: ZERO = 0.0D0, SCALEDIST = 1.2D0
    CHARACTER(LEN=6)    :: J1_STRING

  !1. Bring orientation between first and second loose nucleotide approx. parallel to last base pair
    !First find the vector parallel to the last base pair (approximately)
    BPID = 0
    BPVEC(:) = ZERO
    DO I=1,NUCID-1
       IF (BP_CURR(NUCID-1,I)) THEN
          NUC1AT = LIST_NUCL(I)%LATOM     
          NUC2AT = LIST_NUCL(NUCID-1)%LATOM
          BPID=I
          BPVEC(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
          BPVEC(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
          BPVEC(3) = X(3*NUC2AT) - X(3*NUC1AT)
          EXIT
       ENDIF
    ENDDO
    IF (VEC_LEN(BPVEC).EQ.0.0D0) THEN
       WRITE(MYUNIT, '(A)') " bpmoves> vector for last base pair is zero, no move attempted"
       RETURN
    END IF
    BP_STEPS = BP_STEPS + 1
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_a."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

    !Now find vector between nucid and nucid + 1 (BBVEC1), and vector between nucid - 1 and nucid (BBVEC2)
    BBVEC1(:) = ZERO
    BBVEC2(:) = ZERO
    BBPERP1(:) = ZERO
    BBPERP2(:) = ZERO
    NUC1AT = LIST_NUCL(NUCID)%FATOM     
    NUC2AT = LIST_NUCL(NUCID+1)%FATOM    
    BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
    NUC1AT = LIST_NUCL(NUCID-1)%FATOM     
    NUC2AT = LIST_NUCL(NUCID)%FATOM    
    BBVEC2(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC2(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC2(3) = X(3*NUC2AT) - X(3*NUC1AT)
    !Perpendicular componennt for BBVEC1
    CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
    !Perpendicular component for BPVEC
    CALL PERPVEC(BPVEC,BBVEC2,BBPERP2)
    ! Angle between these two perpendicular components is the rotational angle for the first rotation
    ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
    ROTVEC = ANGLE*VEC_NORM(BBVEC2)
    ! Create rotational matrix
    ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
    ! Select subset of coordinates
    NATTEMP = NATS-LIST_NUCL(NUCID-1)%LATOM !number of atoms in XTEMP
    ALLOCATE(XTEMP(3*NATTEMP))
    FIRSTAT = LIST_NUCL(NUCID)%FATOM
    XTEMP(:) = X(3*FIRSTAT-2:3*NATS)
    ! find new centre (the first atom is of the group is FIRSTAT)
    CENTRE(:) = XTEMP(1:3)
    ! Apply rotation
    DO J=1,NATTEMP
       XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
       XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
       XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
    ENDDO
    ! Apply update of coords to X
    X(3*FIRSTAT-2:3*NATS) = XTEMP(:)
    DEALLOCATE(XTEMP)
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_b."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

  !2. Rotate third free nucleotide to be approximately in line with base pair
    ! define target coordinate for third loose nucleotide
    NUC1AT = LIST_NUCL(NUCID-1)%FATOM
    TARGET2(1) = X(3*NUC1AT-2) + SCALEDIST*BPVEC(1)
    TARGET2(2) = X(3*NUC1AT-1) + SCALEDIST*BPVEC(2)
    TARGET2(3) = X(3*NUC1AT) + SCALEDIST*BPVEC(3)      
    !Now find vector between nucid - 1 and nucid - 2 (BBVEC1), and vector between nucid and nucid -1 (BBVEC2)
    BBVEC1(:) = ZERO
    BBVEC2(:) = ZERO
    BBPERP1(:) = ZERO
    BBPERP2(:) = ZERO
    NUC1AT = LIST_NUCL(NUCID+1)%FATOM     
    NUC2AT = LIST_NUCL(NUCID+2)%FATOM    
    BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
    NUC1AT = LIST_NUCL(NUCID)%FATOM     
    NUC2AT = LIST_NUCL(NUCID+1)%FATOM    
    BBVEC2(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
    BBVEC2(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
    BBVEC2(3) = X(3*NUC2AT) - X(3*NUC1AT)
    ! set vector to TARGET
    NUC1AT = LIST_NUCL(NUCID+1)%FATOM 
    BBTAR(1) = TARGET2(1) - X(3*NUC1AT-2) 
    BBTAR(2) = TARGET2(2) - X(3*NUC1AT-1)
    BBTAR(3) = TARGET2(3) - X(3*NUC1AT)
    !Perpendicular componennt for BBVEC1
    CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
    !Perpendicular component for BBTAR
    CALL PERPVEC(BBTAR,BBVEC2,BBPERP2)
    ! Angle between these two perpendicular components is the rotational angle for the second rotation
    ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
    ROTVEC = ANGLE*VEC_NORM(BBVEC2)
    ! Create rotational matrix
    ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
    ! Select subset of coordinates
    NATTEMP = NATS-LIST_NUCL(NUCID)%LATOM !number of atoms in XTEMP
    ALLOCATE(XTEMP(3*NATTEMP))
    FIRSTAT = LIST_NUCL(NUCID+1)%FATOM
    XTEMP(:) = X(3*FIRSTAT-2:3*NATS)
    ! find new centre
    CENTRE(:) = XTEMP(1:3)
    ! Apply rotation
    DO J=1,NATTEMP
       XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
       XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
       XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
    ENDDO
    ! Apply update of coords to X
    X(3*FIRSTAT-2:3*NATS) = XTEMP(:)
    DEALLOCATE(XTEMP)
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_c."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF

  !3. If we have more than 3 loose nucleotides, rotate the next bond closer to the CoM of rest of the molecule
    IF ((NNUCL-3).GT.NUCID) THEN
       ! first determine CoM of the rest of the molecule
       LASTATOM = LIST_NUCL(NUCID-1)%LATOM
       ALLOCATE(XTEMP(3*LASTATOM))
       XTEMP(:) = X(1:3*LASTATOM)
       COM(:) = ZERO
       DO I=1,LASTATOM
          COM(1) = COM(1) + XTEMP(3*I-2)
          COM(2) = COM(2) + XTEMP(3*I-1)
          COM(3) = COM(3) + XTEMP(3*I)
       ENDDO
       COM(:) = COM(:)/LASTATOM
       DEALLOCATE(XTEMP)
       ! target: vector parallel to middle of last base pair to COM
       NUC1AT = LIST_NUCL(BPID)%LATOM     
       NUC2AT = LIST_NUCL(NUCID-1)%LATOM  
       BPVEC(1) = COM(1) - (X(3*NUC2AT-2) + X(3*NUC1AT-2))/2
       BPVEC(2) = COM(2) - (X(3*NUC2AT-1) + X(3*NUC1AT-1))/2
       BPVEC(3) = COM(3) - (X(3*NUC2AT) + X(3*NUC1AT))/2
       !Now find vector between nucid - 2 and nucid - 3 (BBVEC1), and vector BBVEC1 and BPVEC (BBVEC2)
       BBVEC1(:) = ZERO
       BBVEC2(:) = ZERO
       BBPERP1(:) = ZERO
       BBPERP2(:) = ZERO
       NUC1AT = LIST_NUCL(NUCID+2)%FATOM     
       NUC2AT = LIST_NUCL(NUCID+3)%FATOM    
       BBVEC1(1) = X(3*NUC2AT-2) - X(3*NUC1AT-2)
       BBVEC1(2) = X(3*NUC2AT-1) - X(3*NUC1AT-1)
       BBVEC1(3) = X(3*NUC2AT) - X(3*NUC1AT)
       !rotational reference is between BBVEC1 and TARGET
       BBVEC2(1) = BBVEC1(1) + BPVEC(1)
       BBVEC2(2) = BBVEC1(2) + BPVEC(2)
       BBVEC2(3) = BBVEC1(3) + BPVEC(3)
       !Perpendicular componennt for BBVEC1
       CALL PERPVEC(BBVEC1,BBVEC2,BBPERP1)
       !Perpendicular component for BBTAR
       CALL PERPVEC(BPVEC,BBVEC2,BBPERP2)
       ! Angle between these two perpendicular components is the rotational angle for the first rotation
       ANGLE = VEC_ANGLE(BBPERP1,BBPERP2)
       ROTVEC = ANGLE*VEC_NORM(BBVEC2)
       ! Create rotational matrix
       ROTMAT(:,:)=ROT_AA2MX(ROTVEC(:))
       ! Select subset of coordinates
       NATTEMP = NATS-LIST_NUCL(NUCID+1)%LATOM !number of atoms in XTEMP
       ALLOCATE(XTEMP(3*NATTEMP))
       FIRSTAT = LIST_NUCL(NUCID+2)%FATOM
       XTEMP(:) = X(3*FIRSTAT-2:3*NATS)
       ! find new centre
       CENTRE(:) = XTEMP(1:3)
       ! Apply rotation
       DO J=1,NATTEMP
          XTEMP(3*J-2:3*J) = VEC_DIFF(CENTRE,XTEMP(3*J-2:3*J))     ! shift coordinates
          XTEMP(3*J-2:3*J) = MATMUL(ROTMAT(:,:),XTEMP(3*J-2:3*J))  ! apply rotation
          XTEMP(3*J-2:3*J) = VEC_DIFF(-CENTRE,XTEMP(3*J-2:3*J))    ! shift coordinates back
       ENDDO
       ! Apply update of coords to X
       X(3*FIRSTAT-2:3*NATS) = XTEMP(:)
       DEALLOCATE(XTEMP)
    ENDIF
    !save a pdb file to check steps (comment as needed)
    WRITE(J1_STRING,'(I6)') BP_STEPS
    IF (PRINTPDB.AND.DEBUG) THEN
       CALL DUMP_PDB(3*NATS,X,"bp_hinge_move_d."//TRIM(ADJUSTL(J1_STRING))//".pdb",.TRUE.)
    ENDIF
    RETURN

  END SUBROUTINE MOVE_TTAIL

  ! find perpendicular component of vector relative to a reference
  SUBROUTINE PERPVEC(INPVEC,REFVEC,PERPVCOMP)
    IMPLICIT NONE
    REAL(KIND = REAL64), INTENT(IN) :: INPVEC(3) !vector to be decomposed
    REAL(KIND = REAL64), INTENT(IN) :: REFVEC(3) !reference vector
    REAL(KIND = REAL64), INTENT(OUT) :: PERPVCOMP(3) !perpendicular vector
    REAL(KIND = REAL64) :: DUMMY(3), INPDUMMY(3), REFDUMMY(3)
    
    ! first get normalised vectors
    INPDUMMY = VEC_NORM(INPVEC)
    REFDUMMY = VEC_NORM(REFVEC)
    ! get dot product
    DUMMY = DOT_PRODUCT(REFDUMMY,INPDUMMY)
    ! get perp component
    PERPVCOMP = VEC_NORM(INPDUMMY - DUMMY*INPDUMMY)
    RETURN
  END SUBROUTINE
  
  SUBROUTINE UPDATE_POT(X)
    IMPLICIT NONE
    REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)
    REAL(KIND = REAL64) :: GRADDUMMY(3*NATOMS), EDUMMY
    
    CALL POTENTIAL(X , GRADDUMMY , EDUMMY, .TRUE. , .FALSE.)
    RETURN
  END SUBROUTINE UPDATE_POT
END MODULE BP_MOVES_MOD
