!> @file
!> Contains MOD_NONBONDED to calculate non-bonded contributions

!> Module to calculate all non-bonded contributions including the hydrogen-bonding, excluded volume and stacking energy
MODULE MOD_NONBONDED
   USE PREC_HIRE
   USE NBDEFS

   !> Pair-saturation energy of the last E_NONBONDED call
   REAL(KIND = REAL64) :: ESAT_LAST = 0.0D0
   !> Helix-cooperativity energy of the last E_NONBONDED call
   REAL(KIND = REAL64) :: ECOOP_LAST = 0.0D0
   !> Per-class stacking energy of the last E_NONBONDED call (diagnostic; sums to ESTAK),
   !> used to zero-cost reweight per-class stacking scale factors from logged trajectories
   REAL(KIND = REAL64) :: ESTAK_PYRPUR_LAST = 0.0D0
   REAL(KIND = REAL64) :: ESTAK_PURPUR_LAST = 0.0D0
   REAL(KIND = REAL64) :: ESTAK_PYRPYR_LAST = 0.0D0

   ! Hydrogen-bonded pairs of the current call, kept for the many-body terms
   !> number of stored pairs
   INTEGER, PRIVATE :: NPAIRS = 0
   !> residue indices (i < j) of each pair
   INTEGER, ALLOCATABLE, PRIVATE :: PRES(:,:)
   !> last particle of base i and of base j of each pair
   INTEGER, ALLOCATABLE, PRIVATE :: PATM(:,:)
   !> total and canonical (cWW entry) hydrogen-bond energy of each pair
   REAL(KIND = REAL64), ALLOCATABLE, PRIVATE :: PE(:), PEC(:)
   !> forces of each pair: (:,:,1,p) on base i, (:,:,2,p) on base j, (:,:,3:4,p) the same for the cWW entry
   REAL(KIND = REAL64), ALLOCATABLE, PRIVATE :: PF(:,:,:,:)
   !> pair index of residues (i,j), 0 if no hydrogen bond was stored
   INTEGER, ALLOCATABLE, PRIVATE :: PIDX(:,:)
   !> chain index of each residue (a chain starts at a residue whose first particle is not P)
   INTEGER, ALLOCATABLE, PRIVATE :: RESCHAIN(:)

   CONTAINS
      !> Routine to calaculate non-bonded energy
      !> @brief
      !>
      !> This subroutine calculates the non-bonded terms for all residues.\n
      !> An iteration of all base pairs takes place. Firstly, if the bases are far apart, no further calculations are conducted.\n
      !> Currently, the cutoff value for this is 20 Angstrom and hard-coded.\n
      !> If the bases are close enough, the hydrogen bonding and stacking energies are computed first.\n
      !> Then the excluded volume between particles is calaculated.\n
      !> With the many-body hydrogen-bond terms switched on (KSAT or EPSCOOP non-zero), the hydrogen-bond
      !> forces are collected first and added after the loop by MANYBODY_HB; the saturation and
      !> cooperativity energies are left in ESAT_LAST and ECOOP_LAST.
      !>
      !> @warning The cutoff for calculations is hard coded here to the distance squared being less than 400.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[out] F - gradient array for the non-bonded contributions
      !> @param[out] EHHB - hydrogen bonding energy
      !> @param[out] ESTAK - base stacking energy
      !> @param[out] EVDW - excluded volume interactions
      !>
      !> @see RNA_BB
      !> @see RNA_STACKV2
      !> @see ENERGY_EXCLV
      SUBROUTINE E_NONBONDED(NOPT, X, F, EHHB, ESTAK, EVDW)
         USE UTILS_IO, ONLY: GETUNIT
         USE MOD_HBONDS, ONLY: ENERGY_HB, MINSEP_HB, KSAT, EPSCOOP
         USE MOD_EXCLV, ONLY: ENERGY_EXV
         USE MOD_BASESTACKING, ONLY: NA_STACKV2
         USE NAPARAMS, ONLY: BOCC, BTYPE, BLIST, RCUT2_EXCLV, NBCUT
         USE VAR_DEFS, ONLY: NRES, RESSTART, RESFINAL, RESTYPES, IAC
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force
         REAL(KIND = REAL64), INTENT(OUT) :: EHHB, ESTAK, EVDW

         INTEGER :: I, J           !Iteration variables - residue indices
         INTEGER :: K, L           !Atom indices
         INTEGER :: TYPEI, TYPEJ   !Type of residue (RNA, DNA, protein, ...)
         INTEGER :: TI, TJ         !ID of residue (A, G, ...)
         INTEGER :: TK, TL         !Atom type
         REAL(KIND = REAL64) :: A(3), DA2, DF, DX(3), DCORR
         REAL(KIND = REAL64) :: THIS_EHB, THIS_ESTAK, THIS_EVDW, NB
         LOGICAL :: HBEXIST
         INTEGER :: IA, JA         !last particles of the two bases in a hydrogen bond
         REAL(KIND = REAL64) :: ECWW, FIHB(3,3), FJHB(3,3), FICWW(3,3), FJCWW(3,3)
         LOGICAL :: MANYBODYT

         INTEGER :: STACKUNIT, HBUNIT, EXCLVUNIT

         STACKUNIT = GETUNIT()

#ifdef FOR_ANALYSIS
         OPEN(STACKUNIT, FILE="Dbg_Stacking.dat", STATUS='UNKNOWN')
         HBUNIT = GETUNIT()
         OPEN(HBUNIT, FILE="Dbg_Hbonding_E.dat", STATUS='UNKNOWN')
         EXCLVUNIT = GETUNIT()
         OPEN(EXCLVUNIT, FILE="Dbg_ExcludedV.dat", STATUS='UNKNOWN')
#endif
         EHHB = 0.0D0
         ESTAK = 0.0D0
         EVDW = 0.0D0
         F(1:NOPT) = 0.0D0
         BOCC(1:NRES) = 0

         MANYBODYT = (KSAT.NE.0.0D0).OR.(EPSCOOP.NE.0.0D0)
         ESAT_LAST = 0.0D0
         ECOOP_LAST = 0.0D0
         ESTAK_PYRPUR_LAST = 0.0D0
         ESTAK_PURPUR_LAST = 0.0D0
         ESTAK_PYRPYR_LAST = 0.0D0
         NPAIRS = 0
         IF (.NOT.ALLOCATED(RESCHAIN)) THEN
            CALL SETUP_HB_PAIRS()
         ELSE IF (SIZE(RESCHAIN).NE.NRES) THEN
            CALL SETUP_HB_PAIRS()
         ENDIF

         DO I=1,NRES-1
            DO J=I+1,NRES
               !QUERY: Why are we using first to last here?
               ! Wouldn't it be better to use the distance first to first or last to last?
               !K = BLIST(I)
               !L = BLIST(J-1)+1
               !replaced Blist with RESSTART and RESFINAL to get first and last atom id for all res

               K = RESFINAL(I)
               L = RESSTART(J)
               A(1:3) = X(3*K-2:3*K) - X(3*L-2:3*L)
               DA2 = DOT_PRODUCT(A,A)
               !QUERY: This really should be a variable, not a magic number!
               !Check residues are close enough for interactions
               IF (DA2 .GT. NBCUT) THEN
                  CYCLE
               ENDIF
               !Now calculate interactions between residues
               TYPEI = RESTYPES(I)
               TYPEJ = RESTYPES(J)
               TI = BTYPE(I)
               TJ = BTYPE(J)

               !If restype is 0 it is RNA and if it is 1 it is DNA
               !Currently we can do RNA-RNA or DNA-DNA interactions
               IF (((TYPEI.EQ.0).OR.(TYPEI.EQ.1)).AND.((TYPEJ.EQ.0).OR.(TYPEJ.EQ.1))) THEN
                  !Hydrogen bonding between nucleotides, except between close neighbours of one chain
                  THIS_EHB = 0.0D0
                  IF ((MINSEP_HB.LE.0).OR.(RESCHAIN(I).NE.RESCHAIN(J)).OR.((J-I).GE.MINSEP_HB)) THEN
                     CALL ENERGY_HB(I, J, TYPEI, TYPEJ, NOPT, X, THIS_EHB, HBEXIST, IA, JA, FIHB, FJHB, &
                                    ECWW, FICWW, FJCWW)
                     EHHB = EHHB + THIS_EHB
                     IF (HBEXIST) THEN
                        IF (MANYBODYT) THEN
                           CALL STORE_HB_PAIR(I, J, IA, JA, THIS_EHB, ECWW, FIHB, FJHB, FICWW, FJCWW)
                        ELSE
                           CALL ADD_BASE_FORCES(NOPT, F, IA, JA, FIHB, FJHB, 1.0D0)
                        ENDIF
                     ENDIF
                  ENDIF
                  !Stacking energy
                  !CALL NA_STACKV(NOPT, BLIST(I),BLIST(J),TI,TJ,F,X,THIS_ESTAK,STACKUNIT)  !old Stacking
                  CALL NA_STACKV2(NOPT, BLIST(I),BLIST(J),TI,TJ,F,X,THIS_ESTAK)
                  ESTAK = ESTAK + THIS_ESTAK
                  IF ((TI.LT.3.AND.TJ.GT.2).OR.(TJ.LT.3.AND.TI.GT.2)) THEN
                     ESTAK_PYRPUR_LAST = ESTAK_PYRPUR_LAST + THIS_ESTAK
                  ELSE IF (TI.LT.3.AND.TJ.LT.3) THEN
                     ESTAK_PURPUR_LAST = ESTAK_PURPUR_LAST + THIS_ESTAK
                  ELSE
                     ESTAK_PYRPYR_LAST = ESTAK_PYRPYR_LAST + THIS_ESTAK
                  ENDIF
#ifdef FOR_ANALYSIS
                  IF (ABS(THIS_EHB) .GT. 0.3D0) THEN
                     WRITE(HBUNIT,'(4I6,2F15.7)') I, J, TI, TJ, DSQRT(DA2), THIS_EHB
                  ENDIF
#endif

               !If restype is 2 for both this is protein-protein
               ELSEIF ((TYPEI.EQ.2).AND.(TYPEJ.EQ.2)) THEN
                  CYCLE
               ELSE
                  CYCLE
               ENDIF
               !Now calculate excluded volume interactions
               DO K = RESSTART(I),RESFINAL(I)
                  DO L = RESSTART(J),RESFINAL(J)
                     TK = IAC(K)
                     TL = IAC(L)
                     !make sure we skip the bonded particles for neighbouring res
                     IF ((J-I).EQ.1) THEN
                        IF ((TYPEI.EQ.0).AND.(TYPEJ.EQ.0)) THEN
                           !For RNA ignore O3-P
                           IF ((TK.EQ.1).AND.(TL.EQ.3)) CYCLE
                        ELSEIF ((TYPEI.EQ.1).AND.(TYPEJ.EQ.1)) THEN
                           !For DNA do the same
                           IF ((TK.EQ.1).AND.(TL.EQ.3)) CYCLE
                        ENDIF
                     ENDIF
                     !additionally skip neighbouring nucleotide bb/sugar to bb/sugar distances (types 1,2,3,4)
                     IF ((J-I).EQ.1) THEN
                        IF ((TYPEI.EQ.0).AND.(TYPEJ.EQ.0)) THEN
                           IF ((TK.LE.7).AND.(TL.LE.7)) CYCLE
                        END IF
                     ENDIF
                     A(1:3) = X(3*K-2:3*K) - X(3*L-2:3*L)
                     DA2 = DOT_PRODUCT(A,A)
                     !Skip if the distance is too large
                     IF (DA2.GT.RCUT2_EXCLV) CYCLE
                     DCORR = 1.0D0
                     IF (ABS(I-J).EQ.1) DCORR = 0.5D0
                     CALL ENERGY_EXV(NOPT, X, F, K, L, TK, TL, DA2, DCORR, THIS_EVDW)

                     EVDW = EVDW + THIS_EVDW
#if FOR_ANALYSIS
                     WRITE(EXCLVUNIT,'(4I6,4F15.7)') I, J , K, L, DSQRT(DA2), NBCT2(TK,TL), THIS_EVDW
#endif

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         IF (MANYBODYT.AND.(NPAIRS.GT.0)) CALL MANYBODY_HB(NOPT, F)
#ifdef FOR_ANALYSIS
         CLOSE(HBUNIT)
         CLOSE(STACKUNIT)
         CLOSE(EXCLVUNIT)
#endif
      END SUBROUTINE E_NONBONDED

      !> Many-body hydrogen-bond terms on the pairs stored during E_NONBONDED
      !> @brief
      !>
      !> With the pair energy E_p (<= 0) and the ideal-pair energy EREF of its base types, the occupancy
      !> is o_p = -E_p/EREF and the load of a base n_i = sum over its pairs of o_p.\n
      !> Saturation:   ESAT = KSAT * sum_i max(0, n_i - NSAT0)**2\n
      !> Cooperativity: ECOOP = -EPSCOOP * sum over stacked pair steps (i,j)-(i+1,j-1) of sigma_p*sigma_q,
      !> with sigma = s/(1+s**m)**(1/m) of the cWW occupancy s = -E_cWW/EREF.\n
      !> Both depend on the coordinates only through the pair energies, so the force of each pair is its
      !> hydrogen-bond force times W_p = 1 + dESAT/dE_p, plus its cWW force times C_p = dECOOP/dE_cWW,p.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[inout] F - forces, the scaled hydrogen-bond forces are added
      SUBROUTINE MANYBODY_HB(NOPT, F)
         USE MOD_HBONDS, ONLY: KSAT, NSAT0, EPSCOOP, SIGM, COOPSCALE
         USE VAR_DEFS, ONLY: NRES
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NOPT
         REAL(KIND = REAL64), INTENT(INOUT) :: F(NOPT)

         INTEGER :: P, Q, I, J
         REAL(KIND = REAL64) :: NOCC(NRES), HP(NRES), ER(NPAIRS), SIG(NPAIRS), DSIG(NPAIRS)
         REAL(KIND = REAL64) :: OC, T, W, CF, NEIGH

         NOCC(:) = 0.0D0
         DO P=1,NPAIRS
            I = PRES(1,P)
            J = PRES(2,P)
            ER(P) = PAIR_EREF(I,J)
            SIG(P) = 0.0D0
            DSIG(P) = 0.0D0
            IF (ER(P).LE.0.0D0) CYCLE
            OC = MAX(0.0D0, -PE(P))/ER(P)
            NOCC(I) = NOCC(I) + OC
            NOCC(J) = NOCC(J) + OC
            IF (PEC(P).LT.0.0D0) THEN
               OC = -PEC(P)/ER(P)
               T = 1.0D0 + OC**SIGM
               SIG(P) = OC*T**(-1.0D0/SIGM)
               DSIG(P) = T**(-1.0D0/SIGM - 1.0D0)
            ENDIF
         ENDDO

         HP(:) = 0.0D0
         DO I=1,NRES
            IF (NOCC(I).GT.NSAT0) THEN
               ESAT_LAST = ESAT_LAST + KSAT*(NOCC(I) - NSAT0)**2
               HP(I) = 2.0D0*KSAT*(NOCC(I) - NSAT0)
            ENDIF
         ENDDO

         DO P=1,NPAIRS
            IF (SIG(P).EQ.0.0D0) CYCLE
            Q = HELIX_NEIGHBOUR(P, 1)
            IF (Q.GT.0) ECOOP_LAST = ECOOP_LAST - EPSCOOP*COOPSCALE*SIG(P)*SIG(Q)
         ENDDO

         DO P=1,NPAIRS
            W = 1.0D0
            CF = 0.0D0
            IF (ER(P).GT.0.0D0) THEN
               ! dESAT/dE_p = -(h'(n_i) + h'(n_j))/EREF
               W = 1.0D0 - (HP(PRES(1,P)) + HP(PRES(2,P)))/ER(P)
               IF ((EPSCOOP.NE.0.0D0).AND.(DSIG(P).GT.0.0D0)) THEN
                  NEIGH = 0.0D0
                  Q = HELIX_NEIGHBOUR(P, 1)
                  IF (Q.GT.0) NEIGH = NEIGH + SIG(Q)
                  Q = HELIX_NEIGHBOUR(P, -1)
                  IF (Q.GT.0) NEIGH = NEIGH + SIG(Q)
                  ! dECOOP/dE_cWW,p = EPSCOOP*sigma'(s_p)*(sigma of both neighbouring steps)/EREF
                  CF = EPSCOOP*COOPSCALE*DSIG(P)*NEIGH/ER(P)
               ENDIF
            ENDIF
            CALL ADD_BASE_FORCES(NOPT, F, PATM(1,P), PATM(2,P), PF(:,:,1,P), PF(:,:,2,P), W)
            IF (CF.NE.0.0D0) CALL ADD_BASE_FORCES(NOPT, F, PATM(1,P), PATM(2,P), PF(:,:,3,P), PF(:,:,4,P), CF)
         ENDDO

         DO P=1,NPAIRS
            PIDX(PRES(1,P),PRES(2,P)) = 0
         ENDDO
      END SUBROUTINE MANYBODY_HB

      !> Energy of an ideal pair for residues I and J, 0 unless both are RNA of a known base type
      FUNCTION PAIR_EREF(I, J) RESULT(ER)
         USE MOD_HBONDS, ONLY: EREF_HB
         USE NAPARAMS, ONLY: BTYPE
         USE VAR_DEFS, ONLY: RESTYPES
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: I, J
         REAL(KIND = REAL64) :: ER

         ER = 0.0D0
         IF ((RESTYPES(I).NE.0).OR.(RESTYPES(J).NE.0)) RETURN
         IF ((BTYPE(I).LT.1).OR.(BTYPE(I).GT.4).OR.(BTYPE(J).LT.1).OR.(BTYPE(J).GT.4)) RETURN
         ER = EREF_HB(BTYPE(I),BTYPE(J))
      END FUNCTION PAIR_EREF

      !> Stored pair one helix step away from pair P: DIR = 1 gives (i+1, j-1), DIR = -1 gives (i-1, j+1);
      !> 0 if that step leaves either chain or no hydrogen bond was stored for it
      FUNCTION HELIX_NEIGHBOUR(P, DIR) RESULT(Q)
         USE VAR_DEFS, ONLY: NRES
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: P, DIR
         INTEGER :: Q, I, J

         Q = 0
         I = PRES(1,P) + DIR
         J = PRES(2,P) - DIR
         IF ((I.LT.1).OR.(J.GT.NRES).OR.(I.GE.J)) RETURN
         IF (RESCHAIN(I).NE.RESCHAIN(PRES(1,P))) RETURN
         IF (RESCHAIN(J).NE.RESCHAIN(PRES(2,P))) RETURN
         Q = PIDX(I,J)
      END FUNCTION HELIX_NEIGHBOUR

      !> Store one hydrogen-bonded pair for MANYBODY_HB, growing the arrays when needed
      SUBROUTINE STORE_HB_PAIR(I, J, IA, JA, E, EC, FI, FJ, FIC, FJC)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: I, J, IA, JA
         REAL(KIND = REAL64), INTENT(IN) :: E, EC, FI(3,3), FJ(3,3), FIC(3,3), FJC(3,3)
         INTEGER, ALLOCATABLE :: ITMP(:,:)
         REAL(KIND = REAL64), ALLOCATABLE :: RTMP(:), FTMP(:,:,:,:)
         INTEGER :: NMAX

         NMAX = SIZE(PE)
         IF (NPAIRS.EQ.NMAX) THEN
            ALLOCATE(ITMP(2,2*NMAX))
            ITMP(:,1:NMAX) = PRES
            CALL MOVE_ALLOC(ITMP, PRES)
            ALLOCATE(ITMP(2,2*NMAX))
            ITMP(:,1:NMAX) = PATM
            CALL MOVE_ALLOC(ITMP, PATM)
            ALLOCATE(RTMP(2*NMAX))
            RTMP(1:NMAX) = PE
            CALL MOVE_ALLOC(RTMP, PE)
            ALLOCATE(RTMP(2*NMAX))
            RTMP(1:NMAX) = PEC
            CALL MOVE_ALLOC(RTMP, PEC)
            ALLOCATE(FTMP(3,3,4,2*NMAX))
            FTMP(:,:,:,1:NMAX) = PF
            CALL MOVE_ALLOC(FTMP, PF)
         ENDIF
         NPAIRS = NPAIRS + 1
         PRES(1,NPAIRS) = I
         PRES(2,NPAIRS) = J
         PATM(1,NPAIRS) = IA
         PATM(2,NPAIRS) = JA
         PE(NPAIRS) = E
         PEC(NPAIRS) = EC
         PF(:,:,1,NPAIRS) = FI
         PF(:,:,2,NPAIRS) = FJ
         PF(:,:,3,NPAIRS) = FIC
         PF(:,:,4,NPAIRS) = FJC
         PIDX(I,J) = NPAIRS
      END SUBROUTINE STORE_HB_PAIR

      !> Chain index of every residue and the pair-storage arrays
      SUBROUTINE SETUP_HB_PAIRS()
         USE VAR_DEFS, ONLY: NRES, RESSTART, IGRAPH
         IMPLICIT NONE
         INTEGER :: I, NCH

         IF (ALLOCATED(RESCHAIN)) DEALLOCATE(RESCHAIN, PIDX, PRES, PATM, PE, PEC, PF)
         ALLOCATE(RESCHAIN(NRES), PIDX(NRES,NRES))
         ALLOCATE(PRES(2,4*NRES), PATM(2,4*NRES), PE(4*NRES), PEC(4*NRES), PF(3,3,4,4*NRES))
         PIDX(:,:) = 0
         NCH = 0
         DO I=1,NRES
            IF ((I.EQ.1).OR.(TRIM(ADJUSTL(IGRAPH(RESSTART(I)))).NE.'P')) NCH = NCH + 1
            RESCHAIN(I) = NCH
         ENDDO
      END SUBROUTINE SETUP_HB_PAIRS

      !> Add S times the hydrogen-bond forces of one pair to F (FI acts on IA-2:IA, FJ on JA-2:JA)
      SUBROUTINE ADD_BASE_FORCES(NOPT, F, IA, JA, FI, FJ, S)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NOPT, IA, JA
         REAL(KIND = REAL64), INTENT(INOUT) :: F(NOPT)
         REAL(KIND = REAL64), INTENT(IN) :: FI(3,3), FJ(3,3), S
         INTEGER :: IDX, ID

         DO IDX = 1,3
            ID = IA - IDX + 1
            F((3*ID-2):(3*ID)) = F((3*ID-2):(3*ID)) + S*FI(:,IDX)
            ID = JA - IDX + 1
            F((3*ID-2):(3*ID)) = F((3*ID-2):(3*ID)) + S*FJ(:,IDX)
         ENDDO
      END SUBROUTINE ADD_BASE_FORCES

END MODULE MOD_NONBONDED
