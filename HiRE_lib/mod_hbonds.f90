!> @file 
!> Contains MOD_HBONDS dealing with hydrogen bonding

!> Module contains the functionality to deal with hydrogen bonding contributions to the energy and gradient
MODULE MOD_HBONDS
   USE PREC_HIRE
   USE NBDEFS
   IMPLICIT NONE
   !> Hbond scaling
   REAL(KIND = REAL64) :: EPSHB
   !> titration cutoff for hydrogen bonds
   REAL(KIND = REAL64) :: CTIT
   !> Power for cosine term 
   INTEGER :: P
   !> Gaussian width for length of HB
   REAL(KIND = REAL64) :: Y
   !> Gaussian width for point-plane distance
   REAL(KIND = REAL64) :: GAUSSW
   !> Scaling variable for new planar function
   REAL(KIND = REAL64) :: ALPHA
   !> Interaction scaling
   REAL(KIND = REAL64) :: INTSCALE

   ! Many-body hydrogen-bond terms (applied in MOD_NONBONDED), all off when 0
   !> minimum sequence separation j-i for hydrogen bonds within one chain
   INTEGER :: MINSEP_HB = 0
   !> pair-saturation strength (kcal/mol): KSAT * sum_i max(0, n_i - NSAT0)**2
   REAL(KIND = REAL64) :: KSAT = 0.0D0
   !> pairing occupancy of a base allowed before the saturation penalty starts
   REAL(KIND = REAL64) :: NSAT0 = 1.0D0
   !> helix cooperativity per stacked canonical pair step (kcal/mol): -EPSCOOP * sum sigma_ij sigma_(i+1)(j-1)
   REAL(KIND = REAL64) :: EPSCOOP = 0.0D0
   !> shape exponent m of the capped occupancy sigma(o) = o / (1 + o**m)**(1/m)
   REAL(KIND = REAL64) :: SIGM = 4.0D0
   !> run-time multiplier of EPSCOOP (Hamiltonian replica exchange along the cooperativity, HIRE_INTERFACE:SET_COOP_SCALING)
   REAL(KIND = REAL64) :: COOPSCALE = 1.0D0
   !> table entry of the canonical cWW pair for RNA base types (G 1, A 2, C 3, U 4), 0 for non-canonical type pairs
   INTEGER, PARAMETER :: CWWIDX(4,4) = RESHAPE((/  0,  0, 10, 14, &
                                                   0,  0,  0,  4, &
                                                  10,  0,  0,  0, &
                                                  14,  4,  0,  0 /), (/4,4/))
   !> energy of an ideal pair of RNA base types: EPSHB * s(cWW), or EPSHB * max(s) for non-canonical type pairs
   REAL(KIND = REAL64) :: EREF_HB(4,4)

   CONTAINS

      !> Routine to set module variables based on scale.dat data
      SUBROUTINE SET_HBVARS()
         USE NAPARAMS, ONLY: SCORE_RNA
         USE RNA_HB_PARAMS, ONLY: RS => S, RNPARAM => NPARAM
         INTEGER :: TA, TB
         EPSHB = SCORE_RNA(6)
         CTIT = SCORE_RNA(7)
         P = INT(SCORE_RNA(8))
         Y = SCORE_RNA(9)
         !GAUSSW = SCORE_RNA(73)
         ALPHA = SCORE_RNA(10)
         INTSCALE = SCORE_RNA(11)
         MINSEP_HB = NINT(SCORE_RNA(20))
         KSAT = SCORE_RNA(21)
         NSAT0 = SCORE_RNA(22)
         IF (NSAT0.LE.0.0D0) NSAT0 = 1.0D0
         EPSCOOP = SCORE_RNA(23)
         SIGM = SCORE_RNA(24)
         IF (SIGM.LE.0.0D0) SIGM = 4.0D0
         ! reference energies need the RNA table, filled before this routine (INIT_FROM_MODS)
         DO TA=1,4
            DO TB=1,4
               IF (CWWIDX(TA,TB).GT.0) THEN
                  EREF_HB(TA,TB) = EPSHB*RS(CWWIDX(TA,TB),TA,TB,1,1)
               ELSE
                  EREF_HB(TA,TB) = EPSHB*MAXVAL(RS(1:RNPARAM(TA,TB),TA,TB,1,1))
               ENDIF
            ENDDO
         ENDDO
         IF ((MINSEP_HB.GT.0).OR.(KSAT.NE.0.0D0).OR.(EPSCOOP.NE.0.0D0)) THEN
            WRITE(*,'(A,I3,4(A,F10.4))') " set_hbvars> many-body H-bond terms: MINSEP ", MINSEP_HB, &
               "  KSAT ", KSAT, "  NSAT0 ", NSAT0, "  EPSCOOP ", EPSCOOP, "  SIGM ", SIGM
         ENDIF
      END SUBROUTINE SET_HBVARS

      !> Hydrogen-bond energy of one base pair and its forces on the last three particles of both bases
      !> @brief
      !>
      !> The forces are returned rather than added to F, so that E_NONBONDED can rescale them for the
      !> many-body terms. ECWW, FICWW and FJCWW are the contribution of the canonical cWW table entry alone.
      !>
      !> @param[out] I, JP - last particle of base BI and BJ (the forces act on I-2:I and JP-2:JP)
      SUBROUTINE ENERGY_HB(BI, BJ, MTYPEI, MTYPEJ, NOPT, X, THIS_EHB, HBEXIST, I, JP, FIHB, FJHB, &
                           ECWW, FICWW, FJCWW)
         USE NAPARAMS, ONLY: BTYPE, BP_CURR
         !as the planarityDistEq is zero, we don't need it any longer
!         USE RNA_HB_PARAMS, ONLY: PDEQ_RNA => planarityDistEq
!         USE DNA_HB_PARAMS, ONLY: PDEQ_DNA => planarityDistEq
         USE VAR_DEFS, ONLY: RESFINAL
         USE HB_DEFS, ONLY: SAVE_HB, HBDAT

         INTEGER, INTENT(IN) :: BI, BJ                 ! indices of base I and J
         INTEGER, INTENT(IN) :: MTYPEI, MTYPEJ         ! type of base I and J (RNA or DNA)
         INTEGER, INTENT(IN) :: NOPT                   ! should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    ! input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: THIS_EHB  ! energy contribution
         LOGICAL, INTENT(OUT) :: HBEXIST               ! do we have a hydrogen bond formed?
         INTEGER, INTENT(OUT) :: I, JP                 ! last particle of base I and J
         REAL(KIND = REAL64), DIMENSION(3,3), INTENT(OUT) :: FIHB, FJHB     ! forces on base I and J
         REAL(KIND = REAL64), INTENT(OUT) :: ECWW                           ! energy of the cWW entry
         REAL(KIND = REAL64), DIMENSION(3,3), INTENT(OUT) :: FICWW, FJCWW   ! forces of the cWW entry

         REAL(KIND = REAL64), PARAMETER :: BPTHRESH = 2.3D0 !Energy cutoff for BP in BP_curr

         INTEGER :: J                         ! central atom of base J
         INTEGER :: TI, TJ                    ! base types for I and J - A, C, G, U/T

         INTEGER, PARAMETER :: A = 1 , B = 2
         INTEGER :: IDX, CWWPAR
         REAL(KIND = REAL64) :: ENP1, ENP2, ETEMP, EHHB, REHHB, FTEMP_J(3), FTEMP_I(3)!, DISTEQ
         REAL(KIND = REAL64), DIMENSION(3,3) :: FHB_I, FHB_J, FNP1_I, FNP1_J, FNP2_I, FNP2_J

         ! set variables based on identify of bases
         I = RESFINAL(BI)     ! last atom's index for first base (B1 for A and G, CY for C and U)
         JP = RESFINAL(BJ)    ! last atom's index for base 2
         J = JP - 1            ! central atom's index for base 2
         TI = BTYPE(BI)
         TJ = BTYPE(BJ)
         ! the canonical entry is only tagged for RNA-RNA pairs
         CWWPAR = 0
         IF ((MTYPEI.EQ.0).AND.(MTYPEJ.EQ.0).AND.(TI.GE.1).AND.(TI.LE.4).AND.(TJ.GE.1).AND.(TJ.LE.4)) THEN
            CWWPAR = CWWIDX(TI,TJ)
         ENDIF

         !set forces and energies to zero
         FIHB(:,:) = 0.0D0
         FJHB(:,:) = 0.0D0
         ECWW = 0.0D0
         FICWW(:,:) = 0.0D0
         FJCWW(:,:) = 0.0D0
         THIS_EHB = 0.0D0
         EHHB = 0.0D0
         FHB_I(:,:) = 0.0D0
         FHB_J(:,:) = 0.0D0
         
         ENP1 = 0.0D0
         FNP1_I(:,:) = 0.0D0
         FNP1_J(:,:) = 0.0D0
         ENP2 = 0.0D0
         FNP2_I(:,:) = 0.0D0
         FNP2_J(:,:) = 0.0D0

         ETEMP = 0.0D0
         FTEMP_I(1:3) = 0.0D0
         FTEMP_J(1:3) = 0.0D0
         !QUERY: this is set in the old code - I think there is a problem in how this is all set up and we need to pick either a loop or this without loop
         IDX = 0
         !QUERY: is this the correct one to use here I or J - I guess it is the cross over MTYPEI to TJ
         !determine planarity for I
         !IF (MTYPEI.EQ.0) THEN
         !   DISTEQ = PDEQ_RNA(TJ, 3-IDX)
         !ELSE IF (MTYPEI.EQ.1) THEN
         !   DISTEQ = PDEQ_DNA(TJ, 3-IDX)
         !END IF
         ! Removing PlaneV --> 26/3/2025
         !CALL PlaneV(NOPT, I-B, I-A, I, JP-IDX, X, ENP1, FNP1_I(:,3), FNP1_I(:,2), FNP1_I(:,1), FTEMP_J)!, distEq)
         !determine planarity for J
         !IF (MTYPEJ.EQ.0) THEN
         !   DISTEQ = PDEQ_RNA(TI, 3-IDX)
         !ELSE IF (MTYPEI.EQ.1) THEN
         !   DISTEQ = PDEQ_DNA(TI, 3-IDX)
         !END IF
         ! Removing PlaneV --> 26/3/2025
         !CALL PlaneV(NOPT, JP-B, JP-A, JP, I-IDX, X, ENP2, FNP2_J(:,3), FNP2_J(:,2), FNP1_J(:,1), FTEMP_I)!, distEq)
         !FNP1_J(1:3, IDX + 1) = FNP1_J(1:3, IDX + 1) + FTEMP_J
         !FNP2_I(1:3, IDX + 1) = FNP2_I(1:3, IDX + 1) + FTEMP_I

         CALL HBNEW(BI, BJ, MTYPEI, MTYPEJ, I, TI, J, TJ, NOPT, X, EHHB, HBEXIST, FHB_I, FHB_J, ENP1, ENP2, &
                    CWWPAR, ECWW, FICWW, FJCWW)

         IF (.NOT. HBEXIST) RETURN !at this stage THIS_EHB is still 0.0D0

         !Total energy and regularised energy
         !THIS_EHB = EHHB*(ENP1+ENP2)
         THIS_EHB = EHHB
         REHHB = THIS_EHB/(EPSHB*INTSCALE)
         !update base pairing in BP_CURR
         IF ((ABS(THIS_EHB).GE.BPTHRESH).AND.(ABS(BI-BJ).NE.1)) THEN
            BP_CURR(BI,BJ) = .TRUE.
            BP_CURR(BJ,BI) = .TRUE.
         END IF
         
         ! Additive 
         FIHB = FHB_I !*(ENP1 + ENP2) + EHHB*(FNP1_I + FNP2_I)
         FJHB = FHB_J !*(ENP1 + ENP2) + EHHB*(FNP1_J + FNP2_J)

         !lm759> save Hbond pairs to hbonds.dat
         !IF (SAVE_HB) THEN
         !   IF (ABS(REHHB).GE.1.0D0) THEN
         !      WRITE(HBDAT, '(4i4,4f8.3)') I, JP, BI, BJ, REhhb, Ehhb, Enp1, Enp2
         !   END IF
         !END IF

         ! the forces FIHB/FJHB are added to F by E_NONBONDED (MOD_NONBONDED:ADD_BASE_FORCES)

      END SUBROUTINE ENERGY_HB

      !> Energy and distance between plane and point
      !> @brief
      !>
      !> Computes the distance between one point and the plane defined by 3 other points (i.e. distance(l, plane(i,j,k))).\n
      !> The force and energy contributions are then calculated as well.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] I - index for CG particle i
      !> @param[in] J - index for CG particle j
      !> @param[in] K - index for CG particle k
      !> @param[in] L - index for CG particle l
      !> @param[in] X - coordinates
      !> @param[out] Enewpl - energy for the planarity term
      !> @param[out] FI - forces vector for particle i
      !> @param[out] FJ - forces vector for particle j
      !> @param[out] FK - forces vector for particle k
      !> @param[out] FL - forces vector for particle l                          
      !> @param[in] DEQ - equilibrium distance
      SUBROUTINE PlaneV(NOPT, I, J, K, L, X, Enewpl, FI, FJ, FK, FL)!, DEQ)
         USE VEC_UTILS
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: I,J,K,L                ! particle indices
         INTEGER, INTENT(IN) :: NOPT                   ! number of CG particles
!         REAL(KIND = REAL64), INTENT(IN) :: DEQ     ! equilibrium distance
         
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)   ! coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: Enewpl    ! planarity term
         REAL(KIND = REAL64), INTENT(OUT) :: FI(3), FJ(3), FK(3), FL(3) ! forces on particles

         ! position vectors for particles, and vectors between them
         REAL(KIND = REAL64) :: RI(3), RJ(3), RK(3), RL(3), RIJ(3), RKJ(3), RLJ(3)
         REAL(KIND = REAL64) :: NORMAL(3), DNDQ(3), DIST, DELTA, NNORM, DEDD, DEDD_SCALED, TERM2

         ! define particle positions from coordinates
         RI(1:3) = X(3*I-2:3*I)
         RJ(1:3) = X(3*J-2:3*J)
         RK(1:3) = X(3*K-2:3*K)
         RL(1:3) = X(3*L-2:3*L)
         
         ! get vectors to define plane and distance to point
         RIJ(1:3) = RI(1:3) - RJ(1:3)
         RKJ(1:3) = RK(1:3) - RJ(1:3)
         RLJ(1:3) = RL(1:3) - RJ(1:3)

         ! normal vectors and distance 
         NORMAL = CROSSPRODUCT(RIJ, RKJ)
         NNORM = DSQRT(DOT_PRODUCT(NORMAL, NORMAL))
         DIST = DOT_PRODUCT(NORMAL/NNORM,RLJ)
         
         !IF (DIST.GT.0.0) THEN
         !   DELTA = DIST-DEQ
         !ELSE
         !   DELTA = DIST+DEQ
         !ENDIF
         DELTA = DIST
         ! energy contribution    !OLD: -/3.0 !Add: + !Mul: -
         ENEWPL = INTSCALE * EXP(-(DELTA/GAUSSW)**2)/1.0D0
         ! force 
         DEDD = 2 * ENEWPL * (DELTA/GAUSSW**2)
         DEDD_SCALED = DEDD/NNORM
         TERM2 = DOT_PRODUCT(NORMAL, RLJ)/NNORM**2

         FL(1:3) = DEDD_SCALED*NORMAL(1:3) 
         ! dn / d ri
         FI(1:3) = DEDD_SCALED*(CROSSPRODUCT(RKJ, RLJ)-TERM2*CROSSPRODUCT(RKJ,NORMAL))
         ! dn / d rk
         FK(1:3) = DEDD_SCALED*(CROSSPRODUCT(RLJ, RIJ)-TERM2*CROSSPRODUCT(NORMAL,RIJ))
         ! dn / d rj
         DNDQ = RKJ - RIJ
         FJ(1:3) = DEDD_SCALED*(CROSSPRODUCT(RLJ, DNDQ)-NORMAL-TERM2*CROSSPRODUCT(NORMAL,DNDQ))

      END SUBROUTINE PlaneV

      !> This routine calculates the h-bond energies and forces between two bases.         
      !> @brief
      !>
      !> Hydrogen bonds are calculated based on the orientation od the last three particles in the two bases.\n
      !>  Diagram:\n
      !>      va    ua           ub    vb\n
      !>   a1 -- a2 -- a3 - - b3 -- b2 -- b1\n
      !>              anga   angb\n
      !>
      !> @param[in] BI - index for base i
      !> @param[in] BJ - index for base j
      !> @param[in] IDXA - indices of particle in base A (or I) as reference
      !> @param[in] TYA - type of base A (or I)
      !> @param[in] IDXB - indices of particle in base B (or J) as reference
      !> @param[in] TYB - type of base B (or J)      
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - coordinates
      !> @param[out] EHHB - energy for this hydrogen bond
      !> @param[out] HBEXIST - are the two bases h-bonded?
      !> @param[out] FA - forces array for base A
      !> @param[out] FB - forces array for base b     
      !> @param[in] ENP1 - energy from plane environment to be added for base 1 (or A or I)
      !> @param[in] ENP2 - energy from plane environment to be added for base 2 (or B or J)
      !> @param[in] CWWPAR - table entry of the canonical cWW pair, 0 if none
      !> @param[out] ECWW - energy of the cWW entry alone
      !> @param[out] FACWW, FBCWW - forces of the cWW entry alone
      SUBROUTINE HBNEW(BI, BJ, MTYPEI, MTYPEJ, IDXA, TYA, IDXB, TYB, NOPT, X, EHHB, HBEXIST, &
         FA, FB, ENP1, ENP2, CWWPAR, ECWW, FACWW, FBCWW)
         USE NAPARAMS, ONLY: BPROT, BOCC, RCUT2_HBOND
         USE RNA_HB_PARAMS, ONLY: RCALPAM => CALPAM, RCALPBM => CALPBM, RSALPAM => SALPAM, &
                                  RSALPBM => SALPBM, RDREF => DREF, RS => S, RNPARAM => NPARAM
         USE DNA_HB_PARAMS, ONLY: DCALPAM => CALPAM, DCALPBM => CALPBM, DSALPAM => SALPAM, &
                                  DSALPBM => SALPBM, DDREF => DREF, DS => S, DNPARAM => NPARAM
         USE VEC_UTILS
         IMPLICIT NONE
         
         INTEGER, INTENT(IN) :: BI, BJ                 ! indices of base I (->A) and J (->B)
         !QUERY: A comment in the old routine stated that IDXA and IDXB are the last particles in A and B
         !       But the routine is passed I and J, where J = RESFINAL(B) - 1 
         !       This seems to be accounted for in the definitions of A1 .. B3 later on
         INTEGER, INTENT(IN) :: MTYPEI, MTYPEJ         ! type of base I and J (RNA or DNA)
         INTEGER, INTENT(IN) :: IDXA, IDXB             ! indices of GC particles
         INTEGER, INTENT(IN) :: TYA, TYB               ! types of base A and B
         INTEGER, INTENT(IN) :: NOPT                   ! should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    ! input coordinates
         REAL(KIND = REAL64), INTENT(IN) :: ENP1, ENP2
         REAL(KIND = REAL64), INTENT(OUT) :: EHHB, FA(3,3), FB(3,3)
         LOGICAL, INTENT(OUT) :: HBEXIST
         INTEGER, INTENT(IN) :: CWWPAR
         REAL(KIND = REAL64), INTENT(OUT) :: ECWW, FACWW(3,3), FBCWW(3,3)
         REAL(KIND = REAL64) :: FASAVE(3,3), FBSAVE(3,3)

         REAL(KIND = REAL64), PARAMETER :: REGCUT = 1.0D-3  !Regularisation cutoff  D-7
         
         INTEGER :: QI, QJ        ! charges for I and J, using protonations state
         ! particle positions in A and B
         REAL(KIND = REAL64) :: A1(3), A2(3), A3(3), B1(3), B2(3), B3(3)
         ! vectors between particles in A and B, D... is the norm and ...0 is the normalised vector
         REAL(KIND = REAL64) :: UA(3), UA0(3), DUA, VA(3), DVA
         REAL(KIND = REAL64) :: UB(3), UB0(3), DUB, VB(3), DVB
         REAL(KIND = REAL64) :: RBA(3), DBA, RBA0(3)
         ! crossproducts and their norms and normed crossproducts
         REAL(KIND = REAL64) :: NA(3), NA0(3), DNA, MA0(3), DMA
         REAL(KIND = REAL64) :: NB(3), NB0(3), DNB, MB0(3), DMB
         !further derived variables
         REAL(KIND = REAL64) :: RA(3), DRA, RA0(3), RB(3), DRB, RB0(3)
         REAL(KIND = REAL64) :: COSA, COSB, SINA, SINB
         !local variables to save RNA HB variables from RNA_HB_PARAMS module
         REAL(KIND = REAL64) :: CALPA, SALPA, CALPB, SALPB, SIGHB, STR
         INTEGER :: PAR, NPAR, I !iteration index and limit for iteration
         !variables in the energy and force calculations
         REAL(KIND = REAL64) :: D2, EHHA, VANGL, EHB, dEHB(3)
         REAL(KIND = REAL64) :: ANGA, ANGB, RALPA(3), RALPB(3)
         !new term to account for planarity
         REAL(KIND = REAL64) :: ZA, ZB, FPLAN, DFPLAN, FPLANA(3), FPLANA1(3), FPLANA2(3), FPLANB(3),FPLANB1(3),FPLANB2(3), fanga3(3), fanga2(3), fanga1(3), fangb3(3), fangb2(3), fangb1(3)
         !local replacements for titration globals (TODO: titration needs to be set up properly!) 
         LOGICAL :: use_tit
         INTEGER :: flag_tit 

         use_tit = .FALSE.
         flag_tit = 1

         QI = BPROT(BI)
         QJ = BPROT(BJ)
         EHHB = 0.0D0
         HBEXIST = .FALSE.
         FA(1:3,1:3) = 0.0D0
         FB(1:3,1:3) = 0.0D0
         ECWW = 0.0D0
         FACWW(1:3,1:3) = 0.0D0
         FBCWW(1:3,1:3) = 0.0D0

         !particle positions from coordinates
         A1(1:3) = X((3*IDXA-8):(3*IDXA-6))
         A2(1:3) = X((3*IDXA-5):(3*IDXA-3))
         A3(1:3) = X((3*IDXA-2):(3*IDXA))
         B1(1:3) = X((3*IDXB-5):(3*IDXB-3))
         B2(1:3) = X((3*IDXB-2):(3*IDXB))
         B3(1:3) = X((3*IDXB+1):(3*IDXB+3))
         
         !vectors between particles
         UA(1:3) = A3(1:3) - A2(1:3)
         VA(1:3) = A1(1:3) - A2(1:3)
         UB(1:3) = B3(1:3) - B2(1:3)
         VB(1:3) = B1(1:3) - B2(1:3)     
         RBA(1:3) = A3(1:3) - B3(1:3)

         !check distance between A and B is not too far for HB interactions
         CALL NORMED_VEC(RBA, RBA0, DBA)
         IF (DBA**2 .GE. RCUT2_HBOND) RETURN

         ! get norms and normed vectors for distances
         CALL NORMED_VEC(UA, UA0, DUA)
         CALL NORMED_VEC(UB, UB0, DUB)
         DVA = EUC_NORM(VA)
         DVB = EUC_NORM(VB)
         
         ! get crossproducts
         CALL NORMED_CP2(UA, VA, NA, NA0, DNA)
         CALL NORMED_CP2(UB, VB, NB, NB0, DNB)
         CALL NORMED_CP(NA, UA, MA0, DMA)
         CALL NORMED_CP(NB, UB, MB0, DMB)

         RA(1:3) = -RBA(1:3) - NA0(1:3) * DOT_PRODUCT(-RBA,NA0)
         RB(1:3) = RBA(1:3) - NB0(1:3) * DOT_PRODUCT(RBA,NB0)  
         DRA = EUC_NORM(RA)    
         DRB = EUC_NORM(RB)
         RA0(1:3) = RA(1:3)/DRA
         RB0(1:3) = RB(1:3)/DRB
         
         COSA = DOT_PRODUCT(RA0, UA0)
         SINA = DOT_PRODUCT(RA0, MA0)
         COSB = DOT_PRODUCT(RB0, UB0)
         SINB = DOT_PRODUCT(RB0, MB0)
         
         !add planarity in
         ZA = DOT_PRODUCT(RBA,NA0)/DBA
         ZB = DOT_PRODUCT(RBA,NB0)/DBA
         FPLAN = EXP(-ALPHA*(ZA*ZA+ZB*ZB))

         FPLANA= ZA*(NA0/DBA - RBA*ZA/(DBA**2)  - crossproduct(RBA,VA)/(DBA*DNA) - ZA*(DVA*DVA*UA - dot_product(UA,VA)*VA)/DNA**2) + ZB*(NB0/DBA - RBA*ZB/(DBA**2))

         FPLANA1 =  ZA*((-crossproduct(RBA,UA) + crossproduct(RBA,VA))/(DBA*DNA) + ZA*(DVA*DVA*UA + DUA*DUA*VA - dot_product(UA,VA)*(UA+VA))/DNA**2)

         FPLANA2 = ZA*(+ crossproduct(RBA,UA)/(DBA*DNA) - ZA*(DUA*DUA*VA - dot_product(UA,VA)*UA)/DNA**2)

         FPLANB = ZB*(-NB0/DBA + RBA*ZB/(DBA**2)  - crossproduct(RBA,VB)/(DBA*DNB) - ZB*(DVB*DVB*UB - dot_product(UB,VB)*VB)/DNB**2) - ZA*(NA0/DBA - RBA*ZA/(DBA**2))

         FPLANB1 = ZB*((-crossproduct(RBA,UB) + crossproduct(RBA,VB))/(DBA*DNB) + ZB*(DVB*DVB*UB + DUB*DUB*VB - dot_product(UB,VB)*(UB+VB))/DNB**2)

         FPLANB2 = ZB*(crossproduct(RBA,UB)/(DBA*DNB) - ZB*(DUB*DUB*VB - dot_product(UB,VB)*UB)/DNB**2)


  
         !iteration over all relevant parameters
         IF ((MTYPEI.EQ.0).AND.(MTYPEJ.EQ.0)) THEN
            NPAR = RNPARAM(TYA,TYB)
         ELSE IF ((MTYPEI.EQ.1).AND.(MTYPEJ.EQ.1)) THEN
            NPAR = DNPARAM(TYA,TYB)
         END IF
         DO PAR = 1,NPAR
            !copy HB parameters
            IF ((MTYPEI.EQ.0).AND.(MTYPEJ.EQ.0)) THEN
               SIGHB = RDREF(PAR,TYA,TYB)
               CALPA = RCALPAM(PAR,TYA,TYB)
               SALPA = RSALPAM(PAR,TYA,TYB)
               CALPB = RCALPBM(PAR,TYA,TYB)
               SALPB = RSALPBM(PAR,TYA,TYB)
               STR = RS(PAR,TYA,TYB,QI+1,QJ+1)   !Br2 here is where the WC or non-wc parameters are defined (check)   
            ELSE IF ((MTYPEI.EQ.1).AND.(MTYPEJ.EQ.1)) THEN
               SIGHB = DDREF(PAR,TYA,TYB)
               CALPA = DCALPAM(PAR,TYA,TYB)
               SALPA = DSALPAM(PAR,TYA,TYB)
               CALPB = DCALPBM(PAR,TYA,TYB)
               SALPB = DSALPBM(PAR,TYA,TYB)
               STR = DS(PAR,TYA,TYB,QI+1,QJ+1)   !Br2 here is where the WC or non-wc parameters are defined (check)   
            END IF            

            ! Exponential contribution based on base distance
            D2 = (DBA - SIGHB)/Y
            EHHA = -EPSHB * STR * EXP(-D2**2)

            ! Angular potential contribution (orientation of bases)
            ANGA = COSA*CALPA + SINA*SALPA
            ANGB = COSB*CALPB + SINB*SALPB
            RALPA(1:3) = CALPA*UA0(1:3) + SALPA*MA0(1:3)
            RALPB(1:3) = CALPB*UB0(1:3) + SALPB*MB0(1:3)
            ! warning: P must be an integer otherwise, we cannot get the power for negative ANGA or ANGB
            VANGL = (ANGA*ANGB)**P
            
            !Overall energy for this set of params
            EHB = EHHA*VANGL*FPLAN
            !WRITE(*,*) BI, BJ, STR, IDXA, IDXB, EHHA, EHB
            !d(EHHA)/dX
            dEHB(1:3) = -2.0*EHB*D2/Y*RBA0(1:3)
            

            !first check for size of Ehb
            IF (EHB .GE. REGCUT) CYCLE

            EHHB = EHHB + EHB
            ! the cWW entry's own forces are the change of FA/FB over this iteration
            IF (PAR.EQ.CWWPAR) THEN
               ECWW = EHB
               FASAVE = FA
               FBSAVE = FB
            ENDIF
            !Calculate forces
            
            fa(:,3) = fa(:,3) + fplana2 *2*ALPHA *EHB
            fa(:,2) = fa(:,2) + fplana1 *2*ALPHA *EHB
            fa(:,1) = fa(:,1) + fplana  *2*ALPHA *EHB
            fb(:,1) = fb(:,1) + fplanb  *2*ALPHA *EHB
            fb(:,2) = fb(:,2) + fplanb1 *2*ALPHA *EHB
            fb(:,3) = fb(:,3) + fplanb2 *2*ALPHA *EHB

!             do i = 1, 3
!                if (any(fplana > 100)) then
!                   print*, "Warning: fa(", i, ") fplanea!", bi,bj,fplana
!                endif
!                if (any(fplana1 > 100)) then
!                   print*, "Warning: fa(", i, ") fplanea1!", bi,bj,fplana1
!                endif
!                if (any(fplana2 > 100)) then
!                   print*, "Warning: fa(", i, ") fplanea2!", bi,bj,fplana2
!                endif
!                if (any(fplanb > 100)) then
!                   print*, "Warning: fa(", i, ") fplanea!", bi,bj,fplanb
!                endif
!                if (any(fplanb1 > 100)) then
!                   print*, "Warning: fa(", i, ") fplanea1!", bi,bj,fplanb1
!                endif
!                if (any(fplanb2 > 100)) then
!                endif
!             end do


            fa(:,1) = fa(:,1) - dEhb
            fb(:,1) = fb(:,1) + dEhb

            fanga3 = - Ehb*(p/anga) * (salpa*(crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)))/dma -(crossproduct(ralpa-ra0*anga, ua)*dot_product(-rba, na0)/dna)/dra)
            fanga2 = - Ehb*(p/anga)* (-calpa*(ra0- ua0*cosa)/dua -salpa*( crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)) + &
                     crossproduct(va, crossproduct(ua, ra0-sina*ma0)) + crossproduct(ra0-sina*ma0, na) )/dma -(crossproduct(ua-va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)
            fanga1 = - Ehb*(p/anga)* (calpa*(ra0- ua0*cosa)/dua + salpa*(crossproduct(va, crossproduct(ua, ra0-sina*ma0)) + crossproduct(ra0-sina*ma0, na))/dma &
                     -(ralpa-ra0*anga+crossproduct(va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)- Ehb*(p/angb)*(ralpb-rb0*angb)/drb

            fangb1 = - Ehb*(p/angb)* (calpb*(rb0- ub0*cosb)/dub + salpb*(crossproduct(vb, crossproduct(ub, rb0-sinb*mb0)) + crossproduct(rb0-sinb*mb0, nb))/dmb &
                     -(ralpb-rb0*angb+crossproduct(vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)- Ehb*(p/anga)*(ralpa-ra0*anga)/dra
            fangb2 = - Ehb*(p/angb)* (-calpb*(rb0- ub0*cosb)/dub -salpb*( crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)) +  crossproduct(vb, crossproduct(ub, rb0-sinb*mb0)) + &
                     crossproduct(rb0-sinb*mb0, nb) )/dmb -(crossproduct(ub-vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)
            fangb3 = - Ehb*(p/angb)* (salpb*(crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)))/dmb -(crossproduct(ralpb-rb0*angb, ub)*dot_product(rba, nb0)/dnb)/drb)

            fa(:,3) = fa(:,3) + fanga3
            fa(:,2) = fa(:,2) + fanga2
            fa(:,1) = fa(:,1) + fanga1
            fb(:,1) = fb(:,1) + fangb1
            fb(:,2) = fb(:,2) + fangb2
            fb(:,3) = fb(:,3) + fangb3
            IF (PAR.EQ.CWWPAR) THEN
               FACWW = FA - FASAVE
               FBCWW = FB - FBSAVE
            ENDIF

            
            
!             do i = 1, 3
!                if (any(fanga1 > 100)) then
!                   print*, "Warning: fa(", i, ") fanga1!", bi,bj, fanga1, dra, ZA, ZB, ZA*ZA+ZB*ZB, fplan
!                endif
!                if (any(fanga2 > 100)) then
!                   print*, "Warning: fa(", i, ") fanga2!", bi,bj,fanga2, dra, ZA, ZB, ZA*ZA+ZB*ZB, fplan
!                endif
!                if (any(fanga3 > 100)) then
!                   print*, "Warning: fa(", i, ") fanga3!", bi,bj,fanga3, dra, ZA, ZB, ZA*ZA+ZB*ZB, fplan
!                endif
!                if (any(fangb1 > 100)) then
!                   print*, "Warning: fa(", i, ") fangb1!", bi,bj,fangb1, drb, ZA, ZB, ZA*ZA+ZB*ZB, fplan
!                endif
!                if (any(fangb2 > 100)) then
!                   print*, "Warning: fa(", i, ") fangb2!", bi,bj,fangb2, drb, ZA, ZB, ZA*ZA+ZB*ZB, fplan
!                endif
!                if (any(fangb3 > 100)) then
!                   print*, "Warning: fa(", i, ") fangb3!", bi,bj,fangb3, drb, ZA, ZB,  ZA*ZA+ZB*ZB, fplan
!                endif
!             end do

        ENDDO
        IF (EHHB .GE. REGCUT) RETURN      !! TEST HB EXISTANCE 06-04-2012
         HBEXIST = .TRUE.
      END SUBROUTINE HBNEW

END MODULE MOD_HBONDS
