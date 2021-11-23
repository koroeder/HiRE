MODULE MOD_HBONDS
  USE PREC_HIRE
  USE NBDEFS
  IMPLICIT NONE
  REAL(KIND = REAL64) :: EPSHB, CTIT, P, Y, GAUSSW, INTSCALE  
  
  CONTAINS

   SUBROUTINE SET_HBVARS()
      USE NAPARAMS, ONLY: SCORE_RNA
      EPSHB = SCORE_RNA(23)     
      CTIT = SCORE_RNA(24)
      P = SCORE_RNA(25)
      Y = score_RNA(26)
      GAUSSW = SCORE_RNA(27)
      INTSCALE = SCORE_RNA(28)
   END SUBROUTINE SET_HBVARS 

   SUBROUTINE RNA_BB(BI, BJ, NOPT, X, F, THIS_EHB, HBEXIST)
!-----------------------------------------------------------------------
!>     @details
!>     This function takes care of the Base-Base interaction,
!>     including hydrogen bonding, stacking and cooperativity
!>     the 3 last atoms of each bases are used.
!>
!>     X is the system's coordinates vector
!>     F is the system's force vector
!>     EHHB is the hydrogen bonding energy
!>
!>     I-2 == I-1 == I  - - -  J+1 == J == J-1
!>
!> F : FI3    Fi2   Fi1        Fj1   FJ2   Fj3
!-----------------------------------------------------------------------
      USE NAPARAMS, ONLY: BTYPE, BP_CURR
      USE RNA_HB_PARAMS, ONLY: planarityDistEq
      USE VAR_DEFS, ONLY: RESSTART, RESFINAL
      USE HB_DEFS, ONLY: SAVE_HB, HBDAT
    
      INTEGER, INTENT(IN) :: BI, BJ        ! indices of base I and J
      INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
      REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
      REAL(KIND = REAL64), INTENT(INOUT) :: F(NOPT)   !force from bonds
      REAL(KIND = REAL64), INTENT(OUT) :: THIS_EHB
      LOGICAL, INTENT(OUT) :: HBEXIST      
      
      REAL(KIND = REAL64), PARAMETER :: BPTHRESH = 2.3D0 !Energy cutoff for BP in BP_curr
      
      INTEGER :: I, J                      ! indices for atoms under consideration
      INTEGER :: TI, TJ                    ! base types for I and J

      INTEGER :: IDX, A, B, ID, JP
      REAL(KIND = REAL64) :: ENP1, ENP2, ETEMP, EHHB, REHHB, FTEMP_O(3), DISTEQ
      REAL(KIND = REAL64), DIMENSION(3,3) :: FTEMP, FIHB, FJHB, FIPL, FJPL
      REAL(KIND = REAL64), DIMENSION(3,3) :: FHB_I, FHB_J, FNP1_I, FNP1_J, FNP2_I, FNP2_J
      
      ! set variables based on identify of bases
      I = RESFINAL(BI)     ! last atom's index for first base (B1 for A and G, CY for C and U) 
      JP = RESFINAL(BJ)    ! last atom's index for base 2
      J = JP - 1            ! central atom's index for base 2 
      TI = BTYPE(BI)
      TJ = BTYPE(BJ)

      !set forces and energies to zero
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

      ! back to OLD planarity for indices
      A = 1 !NEW: 3 /OLD: 1 
      B = 2 !NEW: 5 /OLD: 2
!      IF (I < 6) B=4 !NEW planarity

      DO IDX = 0,0 !NEW: 0,0 /OLD: 0,2
         ETEMP = 0.0D0
         FTEMP(:,:) = 0.0D0
         FTEMP_o(:) = 0.0D0
         DISTEQ = planarityDistEq(TJ, 3-IDX)
         CALL RNA_NewPlanev(NOPT, I-b, I-a, I, JP-IDX, X, Etemp, Ftemp(:,3), Ftemp(:,2), Ftemp(:,1), Ftemp_o, distEq)
         Fnp1_i = Fnp1_i + Ftemp
         Fnp1_j(:,idx+1) = Fnp1_j(:,idx+1) + Ftemp_o
         Enp1 = Enp1 + Etemp

         ETEMP = 0.0D0
         FTEMP(:,:) = 0.0D0
         FTEMP_o(:) = 0.0D0
         DISTEQ = planarityDistEq(TI, 3-idx)
         CALL RNA_NewPlanev(NOPT, JP-b, JP-a, JP, I-IDX, X, Etemp, Ftemp(:,3), Ftemp(:,2), Ftemp(:,1), Ftemp_o, distEq)
         Fnp2_j = Fnp2_j + Ftemp
         Fnp2_i(:,idx+1) = Fnp2_i(:,idx+1) + Ftemp_o
         Enp2 = Enp2 + Etemp
      END DO

      !get HB contribution
      CALL RNA_HBNEW(BI, BJ, I, TI, J, TJ, NOPT, X, EHHB, HBEXIST, &
                     FHB_I, FHB_J, ENP1, ENP2)

      IF (.NOT. HBEXIST) RETURN !at this stage THIS_EHB is still 0.0D0

      !Total energy and regularised energy
      !QUERY: which one are we now using - Mult or add? What's the difference?
      THIS_EHB = EHHB*(ENP1+ENP2)         !Mult: * !Add: +
      REHHB = THIS_EHB/(EPSHB*INTSCALE)
      !update base pairing in BP_CURR
      IF ((ABS(THIS_EHB).GE.BPTHRESH).AND.(ABS(BI-BJ).NE.1)) THEN
          BP_CURR(BI,BJ) = .TRUE.
          BP_CURR(BJ,BI) = .TRUE.
      END IF
      
      !
      ! Multiplicative Energy
!      Fihb = Fhb_i*Enp1*Enp2
!      Fjhb = Fhb_j*Enp1*Enp2
!      Fipl = Ehhb*Fnp1_i*Enp2 + Ehhb*Enp1*Fnp2_i
!      Fjpl = Ehhb*Fnp1_j*Enp2 + Ehhb*Enp1*Fnp2_j
      ! Additive Energy
      FIHB = FHB_I*(ENP1 + ENP2)
      FJHB = FHB_J*(ENP1 + ENP2)
      FIPL = EHHB*(FNP1_I + FNP2_I)
      FJPL = EHHB*(FNP1_J + FNP2_J)

      ! OLD planarity !NEW Additive
      FIHB = FIHB + FIPL
      FJHB = FJHB + FJPL

      !lm759> save Hbond pairs to hbonds.dat
      IF (SAVE_HB) THEN
         IF (ABS(REHHB).GE.1.0D0) THEN
            WRITE(HBDAT, '(4i4,4f8.3)') I, JP, BI, BJ, REhhb, Ehhb, Enp1, Enp2
         END IF
      END IF

      ! NEW planarity Multiplicative
!      id = I
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,1) +  Fipl(:,1)
!      id = J +1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,1) +  Fjpl(:,1)
!      id = I-1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,2) 
!      id = I-2
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fihb(:,3)
!      id = J
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,2)
!      id = J-1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjhb(:,3)
!      id = I - a
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fipl(:,2)
!      id = J - a+1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjpl(:,2)
!      id = I - b
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fipl(:,3)
!      id = J - b+1
!      F(id*3-2:id*3)=F(id*3-2:id*3) + Fjpl(:,3)

      ! OLD planarity !NEW Additive
      DO IDX = 1,3
         ID = I - IDX + 1
         F((3*ID-2):(3*ID)) = F((3*ID-2):(3*ID)) + FIHB(:,IDX)
         ID = JP - IDX + 1
         F((3*ID-2):(3*ID)) = F((3*ID-2):(3*ID)) + FJHB(:,IDX)      
      ENDDO
   END SUBROUTINE RNA_BB


   SUBROUTINE RNA_HBNEW(BI, BJ, IDXA, TYA, IDXB, TYB, NOPT, X, EHHB, HBEXIST, &
                        FA, FB, ENP1, ENP2)
!-----------------------------------------------------------------------
!>  @brief
!>  This routine calculates the h-bond energies and forces between two bases.
!
!>  @details
!>  hbexist is a boolean whose value will depend on the presence of an h-bond
!>
!>    Diagram:
!>
!>      va    ua           ub    vb
!>   a1 -- a2 -- a3 - - b3 -- b2 -- b1
!>              anga   angb
!>
!-----------------------------------------------------------------------
      USE NAPARAMS, ONLY: BPROT, BOCC, RCUT2_HB_MCMC_OUT
      USE RNA_HB_PARAMS
      USE VEC_UTILS
      USE HB_DEFS, ONLY :SAVE_HB
      IMPLICIT NONE
      
      INTEGER, INTENT(IN) :: BI, BJ                 ! indices of base I (->A) and J (->B)
      !QUERY: A comment in the old routine stated that IDXA and IDXB are the last particles in A and B
      !       But the routine is passed I and J, where J = RESFINAL(B) - 1 
      !       This seems to be accounted for in the definitions of A1 .. B3 later on
      INTEGER, INTENT(IN) :: IDXA, IDXB             ! indices of GC particles
      INTEGER, INTENT(IN) :: TYA, TYB               ! types of base A and B
      INTEGER, INTENT(IN) :: NOPT                   ! should be 3*NATOMS
      REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    ! input coordinates
      REAL(KIND = REAL64), INTENT(IN) :: ENP1, ENP2
      REAL(KIND = REAL64), INTENT(OUT) :: EHHB, FA(3,3), FB(3,3) 
      LOGICAL, INTENT(OUT) :: HBEXIST
      
      REAL(KIND = REAL64), PARAMETER :: REGCUT = 1.0D-7  !Regularisation cutoff
      
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
      INTEGER :: PAR  !iteration index
      !variables in the energy and force calculations
      REAL(KIND = REAL64) :: D2, EHHA, VANGL, EHB, dEHB(3)
      REAL(KIND = REAL64) :: ANGA, ANGB, RALPA(3), RALPB(3)
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
      IF (DBA**2 .GE. RCUT2_HB_MCMC_OUT) RETURN
      
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
      
      !iteration over all relevant parameters
      DO PAR = 1,NPARAM(TYA,TYB)
         !copy HB parameters
         SIGHB = dREF(PAR,TYA,TYB)
         CALPA = CALPAM(PAR,TYA,TYB)
         SALPA = SALPAM(PAR,TYA,TYB)
         CALPB = CALPBM(PAR,TYA,TYB)
         SALPB = SALPBM(PAR,TYA,TYB)
         STR = S(PAR,TYA,TYB,QI+1,QJ+1)   !Br2 here is where the WC or non-wc parameters are defined (check)

         ! Exponential contribution based on base distance
         D2 = (DBA - SIGHB)/Y
         EHHA = -EPSHB * STR * EXP(-D2**2)
         
         ! Angular potential contribution (orientation of bases)
         ANGA = COSA*CALPA + SINA*SALPA
         ANGB = COSB*CALPB + SINB*SALPB
         RALPA(1:3) = CALPA*UA0(1:3) + SALPA*MA0(1:3)
         RALPB(1:3) = CALPB*UB0(1:3) + SALPB*MB0(1:3)
         VANGL = (ANGA*ANGB)**P
         
         !Overall energy for this set of params
         EHB = EHHA*VANGL
         dEHB(1:3) = -2.0*EHB*D2/Y*RBA0(1:3)
         
         !first check for size of Ehb
         IF (EHB .GE. REGCUT) CYCLE
         !QUERY: the potential document mentions two regularisations (p. 9, eq. 46,47),
         !       one at 1.0D-8 and one at 1.0D-7
         !       I can only find two checks at  1.0D-7
         !       This value is now saved as REGCUT, but should there be another 
         !       check or a different value used?
         
         !TODO: work on titration
         !WARNING: currently use_tit and flag_tit are set as local variables
         !         and the following code is always skipped 
         
         !determine what bases are not free for protonation due to HB          ! Br3 includes planarity
         !QUERY isn't this the wrong place to test this? isn't it EHHB outside this DO loop we want?
         if (use_tit .and. flag_tit .eq. 1) then
            if (abs(EHB*Enp1*Enp2) .ge. ctit) then                             
               ! Br2 need to find a good cutoff, 2.0 is too big --> 1. need to protect more HB
               if (s(par,tya,tyb,1,qj+1) .ne. s(par,tya,tyb,2,qj+1)) then
                  bocc(bi)=1
               endif
               if (s(par,tya,tyb,qi+1,1) .ne. s(par,tya,tyb,qi+1,2)) then
                  bocc(bj)=1
               endif
            endif
         endif
         !end of titration section
         !Update energy
         EHHB = EHHB + EHB   
         !Calculate forces
         fa(:,3) = fa(:,3) - Ehb*(p/anga) * (&
!                   d anga / d ma0
          salpa*(crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)))/dma &
!                   d anga / d ra0
          -(crossproduct(ralpa-ra0*anga, ua)*dot_product(-rba, na0)/dna)/dra)
         fa(:,2) = fa(:,2)        - Ehb*(p/anga)* (&
!                   d anga / d ua0
          -calpa*(ra0- ua0*cosa)/dua &
!                   d anga / d ma0
          -salpa*( crossproduct(ua0, crossproduct(ra0-sina*ma0, ua0)) + &
          crossproduct(va, crossproduct(ua, ra0-sina*ma0)) + & 
          crossproduct(ra0-sina*ma0, na) )/dma &
!                   d anga / d ra0
          -(crossproduct(ua-va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)
         fa(:,1) = fa(:,1) - dEhb - Ehb*(p/anga)* (&
!                   d anga / d ua0
          calpa*(ra0- ua0*cosa)/dua +&
!                   d anga / d ma0
          salpa*(crossproduct(va, crossproduct(ua, ra0-sina*ma0)) + &
          crossproduct(ra0-sina*ma0, na))/dma &
!                   d anga / d ra0
          -(ralpa-ra0*anga+crossproduct(va, ralpa-ra0*anga)*dot_product(-rba, na0)/dna)/dra)&
!                   d angb / d rb0
          - Ehb*(p/angb)*(ralpb-rb0*angb)/drb
         fb(:,1) = fb(:,1) + dEhb - Ehb*(p/angb)* (&
!                   d angb / d ub0
          calpb*(rb0- ub0*cosb)/dub +&
!                   d angb / d mb0
          salpb*(crossproduct(vb, crossproduct(ub, rb0-sinb*mb0)) + &
          crossproduct(rb0-sinb*mb0, nb))/dmb &
!                   d angb / d rb0
          -(ralpb-rb0*angb+crossproduct(vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)&
!                   d anga / d ra0
          - Ehb*(p/anga)*(ralpa-ra0*anga)/dra
         fb(:,2) = fb(:,2)        - Ehb*(p/angb)* (&
!                   d angb / d ub0
          -calpb*(rb0- ub0*cosb)/dub &
!                   d angb / d mb0
          -salpb*( crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)) + &
          crossproduct(vb, crossproduct(ub, rb0-sinb*mb0)) + &
          crossproduct(rb0-sinb*mb0, nb) )/dmb &
!                   d angb / d rb0
          -(crossproduct(ub-vb, ralpb-rb0*angb)*dot_product(rba, nb0)/dnb)/drb)
         fb(:,3) = fb(:,3)        - Ehb*(p/angb)* (&
!                   d angb / d mb0
          salpb*(crossproduct(ub0, crossproduct(rb0-sinb*mb0, ub0)))/dmb &
!                   d angb / d rb0
          -(crossproduct(ralpb-rb0*angb, ub)*dot_product(rba, nb0)/dnb)/drb)

      ENDDO
      IF (EHHB .GE. REGCUT) RETURN      !! TEST HB EXISTANCE 06-04-2012
      HBEXIST = .TRUE.
   END SUBROUTINE RNA_HBNEW

!-------------------------------------------------------------------------
!>    @brief
!>     Computes the distance between one point and the plane defined by 3 other points.
!>     distance(l, plane(i,j,k))
!>     The force and energy contributions are then calculated as well.
!-------------------------------------------------------------------------
   SUBROUTINE RNA_NewPlanev(NOPT, I, J, K, L, X, Enewpl, FI, FJ, FK, FL, distEq)
      USE VEC_UTILS
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: I,J,K,L                ! particle indices
      INTEGER, INTENT(IN) :: NOPT                   ! number of CG particles
      REAL(KIND = REAL64), INTENT(IN) :: DISTEQ     ! equilibrium distance
      
      REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)   ! coordinates
      REAL(KIND = REAL64), INTENT(OUT) :: Enewpl    ! planarity term
      REAL(KIND = REAL64), INTENT(OUT) :: FI(3), FJ(3), FK(3), FL(3) ! forces on particles

      ! position vectors for particles, and vectors between them
      REAL(KIND = REAL64) :: RI(3), RJ(3), RK(3), RL(3), RIJ(3), RKJ(3), RLJ(3)
      REAL(KIND = REAL64) :: NORMAL(3), DNDQ(3), DIST, DELTA, NNORM, DEDD

      
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
      
      ! QUERY: is this correct? Or do we want the absolute value?
      IF (DIST.GT.0.0) THEN
         DELTA = DIST-DISTEQ
      ELSE
         DELTA = DIST+DISTEQ
      ENDIF
      
      ! energy contribution    !OLD: -/3.0 !Add: + !Mul: -
      ENEWPL = INTSCALE * EXP(-(DELTA/GAUSSW)**2)/1.0D0
      ! force 
      DEDD = 2 * ENEWPL * (DELTA/GAUSSW**2)
      FL(1:3) = DEDD*NORMAL(1:3)/NNORM   
      
      ! dn / d ri
      FI(1:3) = dedd*(crossproduct(RKJ, RLJ) - dot_product(normal, RLJ) &
                      * crossproduct(RKJ,normal)/nnorm**2)/nnorm
      ! dn / d rk
      FK(1:3) = dedd*(crossproduct(RLJ, RIJ) - dot_product(normal, RLJ) &
                      * crossproduct(normal,RIJ)/nnorm**2)/nnorm     
      ! dn / d rj
      DNDQ = RKJ - RIJ
      FJ(1:3) = dedd*(crossproduct(RLJ, DNDQ)-normal- dot_product(normal, RLJ) & 
                      * crossproduct(normal,DNDQ)/nnorm**2)/nnorm

   END SUBROUTINE RNA_NewPlanev
   
END MODULE MOD_HBONDS
