MODULE MOD_DIHEDRALS
  USE PREC_HIRE
  INTEGER :: NDIHS    !Number of dihedrals
  INTEGER :: NPTRA    !Number of unique dihedral types

  !Dihedrals: E = PK * (1 + COS(PN * P - PHASE)) with P the current dihedral angle
  REAL(KIND = REAL64), ALLOCATABLE :: PK(:)    !Torsional force constant
  REAL(KIND = REAL64), ALLOCATABLE :: PN(:)    !Number of torsional equilibrium states 
  INTEGER, ALLOCATABLE :: IPN(:)               !Interger value of PN
  REAL(KIND = REAL64), ALLOCATABLE :: PHASE(:) !Torsional phase factor
  REAL(KIND = REAL64), ALLOCATABLE :: GAMC(:)  !Reference values for torsions (cosine)
  REAL(KIND = REAL64), ALLOCATABLE :: GAMS(:)  !Reference values for torsions (sine)
  REAL(KIND = REAL64), ALLOCATABLE :: FMN(:)   !Multiplication used in determining forces for multiple equilibrium angles

  !book-keeping for torsionals
  INTEGER, ALLOCATABLE :: IP(:)                 !Dihedral atom1
  INTEGER, ALLOCATABLE :: JP(:)                 !Dihedral atom2
  INTEGER, ALLOCATABLE :: KP(:)                 !Dihedral atom3
  INTEGER, ALLOCATABLE :: LP(:)                 !Dihedral atom4
  INTEGER, ALLOCATABLE :: ICP(:)                !Type of dihedral

  INTEGER, ALLOCATABLE :: NPHITYPE(:)   !Torsion scaling for backbone (0) or base (1)

  CONTAINS

   SUBROUTINE ALLOC_DIHS()
      CALL DEALLOC_DIHS()
      ALLOCATE(PK(NPTRA), PN(NPTRA), IPN(NPTRA), PHASE(NPTRA), &
               IP(NDIHS), JP(NDIHS), KP(NDIHS), LP(NDIHS), ICP(NDIHS), &
               GAMC(NPTRA), GAMS(NPTRA), FMN(NPTRA), NPHITYPE(NDIHS))
   END SUBROUTINE ALLOC_DIHS 

   SUBROUTINE INIT_DIHPAR()
      USE NUM_DEFS, ONLY: PI
      IMPLICIT NONE
      REAL(KIND = REAL64) :: DUM, DUMS, DUMC
      REAL(KIND = REAL64), PARAMETER :: EPS1 = 1.0D-3
      REAL(KIND = REAL64), PARAMETER :: EPS2 = 1.0D-6
      INTEGER :: I

      DO I = 1,NPTRA
         DUM = PHASE(I)
         IF (DABS(DUM-PI) .LE. EPS1) DUM = SIGN(PI,DUM)
         DUMC = DCOS(DUM)
         DUMS = DSIN(DUM)
         IF(DABS(DUMC) .LE. EPS2) DUMC = 0.0d0
         IF(DABS(DUMS) .LE. EPS2) DUMS = 0.0d0

         GAMC(I) = DUMC*PK(I)
         GAMS(I) = DUMS*PK(I)
 
         FMN(I) = 1.0d0
         IF(PN(I) .LE. 0.0d0) FMN(I) = 0.0d0  
         PN(I) = DABS(PN(I))
         IPN(I) = INT(PN(I)+EPS1)      
      ENDDO
   END SUBROUTINE INIT_DIHPAR

   SUBROUTINE ASSIGN_PHITYPE()
      USE VAR_DEFS, ONLY: IAC   
      IMPLICIT NONE
      INTEGER :: JN, I3, J3, K3T, L3T, IT0, IT1, IT2, IT3
      
      DO JN=1,NDIHS
         I3 = IP(JN)/3 + 1
         J3 = JP(JN)/3 + 1
         K3T = KP(JN)/3 + 1
         L3T = LP(JN)/3 + 1
         IT0 = IAC(I3)
         IT1 = IAC(J3)
         IT2 = IAC(K3T)
         IT3 = IAC(L3T)
         !QUERY: should there be a default value here or should we raise an error?
         IF (IT0.EQ.1) THEN
            !IT1=4 ; IT2=5 ; IT3=(6,8,10,11)
            IF (IT2.EQ.5) THEN
               NPHITYPE(JN) = 2
            !IT1=4 ; IT2=3 ; IT3=2  
            ELSE IF (IT2.EQ.3) THEN
               NPHITYPE(JN) = 4
            ENDIF
         ELSE IF (IT0.EQ.2) THEN
            !IT1=1 ; IT2=4 ; IT3=3
            IF (IT3.EQ.3) THEN
               NPHITYPE(JN) = 6
            !IT1=1 ; IT2=4 ; IT3=5  
            ELSE IF (IT3.EQ.5) THEN
               NPHITYPE(JN) = 7
            ENDIF        
         ELSE IF (IT0.EQ.3) THEN
            !IT1=4 ; IT2=5 ; IT3=(6,8,10,11)
            IF (IT1.EQ.4) THEN
               NPHITYPE(JN) = 3
            !IT1=2 ; IT2=1 ; IT3=4  
            ELSE IF (IT1.EQ.2) THEN
               NPHITYPE(JN) = 8
            ENDIF         
         ELSE IF (IT0.EQ.4) THEN
            !IT1=5 ; IT2=(6,8) ; IT3=(7,9)
            IF (IT1.EQ.5) THEN
               NPHITYPE(JN) = 0
            !IT1=(6,8) ; IT2=(7,9) ; IT3=5  
            ELSE IF ((IT1.EQ.6).OR.(IT1.EQ.8)) THEN
               NPHITYPE(JN) = 1
            !IT1=3 ; IT2=2 ; IT3=1  
            ELSE IF (IT1.EQ.3) THEN
               NPHITYPE(JN) = 9               
            ENDIF         
         ELSE IF (IT0.EQ.5) THEN
            !IT1=4 ; IT2=3 ; IT3=2
            IF (IT1.EQ.4) THEN
               NPHITYPE(JN) = 5         
            ENDIF
         ENDIF
      ENDDO
   END SUBROUTINE ASSIGN_PHITYPE

   SUBROUTINE GET_DIHEDRALS(DIHINFO)
      INTEGER, INTENT(OUT) :: DIHINFO(NDIHS,5)
      INTEGER JN

      DO JN=1,NDIHS
         DIHINFO(JN,1) = IP(JN)/3 + 1
         DIHINFO(JN,2) = JP(JN)/3 + 1
         DIHINFO(JN,3) = IABS(KP(JN))/3 + 1
         DIHINFO(JN,4) = IABS(LP(JN))/3 + 1
         DIHINFO(JN,5) = NPHITYPE(JN)         
      ENDDO
   END SUBROUTINE GET_DIHEDRALS

   SUBROUTINE ENERGY_DIHS(NOPT, X, F, ETORS)
      USE VEC_UTILS
      USE NUM_DEFS, ONLY: PI
      USE NAPARAMS, ONLY: SCORE_RNA
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
      REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
      REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
      REAL(KIND = REAL64), INTENT(OUT) :: ETORS
      !regularisation parameters     
      REAL(KIND = REAL64), PARAMETER :: EPS_R = 1.0D-9
      REAL(KIND = REAL64), PARAMETER :: TM06 = 1.0d-06
      REAL(KIND = REAL64), PARAMETER :: EPS_Z = 1.0d-03
      REAL(KIND = REAL64), PARAMETER :: TM24 = 1.0d-18
      !integers containing atom indices
      INTEGER :: I3, J3, K3, K3T, L3, L3T
      !vectors
      REAL(KIND = REAL64) :: RIJ(3), RKJ(3), RKL(3), RD(3), RG(3)
      REAL(KIND = REAL64) :: rFI(3), rFJ(3), rFK(3), rFL(3) 
      REAL(KIND = REAL64) :: rDC(3), rDC2(3), rDR1(3), rDR2(3), rDR(3)

      INTEGER :: JN, IC, INC
      !TODO: work out the meaning of the variables and comment as much as possible
      !numerous variables
      REAL(KIND = REAL64) :: COSNPs, SINNPs, EXPRB, vEPW, DF0, DF1 
      REAL(KIND = REAL64) :: DFLIM, DUMS, vDF
      REAL(KIND = REAL64) :: LEND, LENG, DOTDG, Z10, Z20, Z12, vFMUL
      REAL(KIND = REAL64) :: CT0, CT1, AP1, vCPHI, vSPHI, COSNP, SINNP
      REAL(KIND = REAL64) :: GMUL(10)
      
      GMUL=(/0.0d+00, 2.0d+00, 0.0d+00, 4.0d+00, 0.0d+00, 6.0d+00, &  
             0.0d+00, 8.0d+00, 0.0d+00, 10.0d+00/) 
      ETORS = 0.0D0
      F(1:NOPT) = 0.0D0

      DO JN = 1,NDIHS
         !assign indices for atoms in dihedral
         I3 = IP(JN)
         J3 = JP(JN)
         K3T = KP(JN)
         L3T = LP(JN)
         K3 = IABS(K3T)
         L3 = IABS(L3T)
         
         !calculate vectors between atoms 
         RIJ = X(I3+1:I3+3) - X(J3+1:J3+3)
         RKJ = X(K3+1:K3+3) - X(J3+1:J3+3)
         RKL = X(K3+1:K3+3) - X(L3+1:L3+3)

         RD = crossproduct(RIJ, RKJ)
         RG = crossproduct(RKL, RKJ)

         !Regularisation from here
         !QUERY: This is described around eqn (10), but the equations are LEN + EPSR**2
         !       Which version is correct?
         LEND = DSQRT(DOT_PRODUCT(RD, RD)+EPS_R**2)
         LENG = DSQRT(DOT_PRODUCT(RG, RG)+EPS_R**2) 

         DOTDG = dot_product(RD,RG)

         Z10 = 1.0d0/LEND
         Z20 = 1.0d0/LENG
         
         IF (EPS_Z .GT. LEND) Z10 = 0.0D0
         IF (EPS_Z .GT. LENG) Z20 = 0.0D0
         Z12 = Z10*Z20
         
         vFMUL = 0.0D0
         IF (Z12 .NE. 0.0d0) vFMUL = 1.0d0

         CT0 = MIN(1.0d0,DOTDG*Z12)
         CT1 = MAX(-1.0d0,CT0)

         AP1 = PI-DSIGN(DACOS(CT1),dot_product(rKJ,crossproduct(rG,rD)))
         vCPHI = -CT1
         vSPHI = DSIN(AP1)

         ! regularised energy is given by:
         ! E = kd*[1 + cos(m*phi_reg)*c0eps + Deps*sin(m*phi_reg)*s0eps] 
         ! This reduces to E = S3*kd*[1 + cos(m*phi - phi0)]
         ! S3 - weight coefficient
         ! Phi0 - equilibrium angle
         ! kd - coupling coefficient
         ! m - number of equilibrium states
         ! These quantities are computed as follows:
         ! phi_reg is pi - sign(dot(rkj,(nj,nk)))*arccos(nj*nk), where nj and nk 
         ! are the normal vectors to the planes ijk and jkl, this variable is AP1
         ! s0eps is sin(phi_reg0)
         
         ! ----- ENERGY AND THE DERIVATIVES WITH RESPECT TO COSPHI -----
         IC = ICP(JN)
         INC = IPN(IC)
         CT0 = PN(IC)*AP1
         COSNP = DCOS(CT0)
         SINNP = DSIN(CT0)

       if(PN(IC).eq.12)then
            !GAMCs=GAMC(IC)/PK(IC)
            !GAMSs=GAMS(IC)/PK(IC)
            COSNPs=DCOS(AP1)
            SINNPs=DSIN(AP1)
            EXPRB=ACOS(GAMC(IC)/PK(IC))*180.d0/PI

            vEPW=(PK(IC)*COSNPs**int(EXPRB))*vFMUL

            if(EXPRB.eq.0)then
               DF0=0.d0
               df1=0.d0
            else
               DF0=PK(IC)*EXPRB*SINNPs*COSNPs**int(EXPRB-1.d0)
               !print*,DF0
               DUMS = vSPHI+SIGN(TM24,vSPHI)
               DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
               
               df1 = df0/dums
               if(tm06.gt.abs(dums)) df1 = dflim
            endif
            !DF0 = -PN(IC)*PK(IC)*(GAMCs*SINNP-GAMSs*COSNP)*(COSNP*GAMCs+SINNP*GAMSs)**(PN(IC)-1.d0)
         else
            vEPW= (PK(IC)+COSNP*GAMC(IC)+SINNP*GAMS(IC))*vFMUL !! might be revised
            DF0 = PN(IC)*(GAMC(IC)*SINNP-GAMS(IC)*COSNP)
            DUMS = vSPHI+SIGN(TM24,vSPHI)
            DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
            
            df1 = df0/dums
            if(tm06.gt.abs(dums)) df1 = dflim
        endif
       
    
        vDF = DF1*vFMUL

        vEPW = vEPW*score_RNA(10)*score_RNA(29+nphitype(JN))   ! score(16) backbone k, score(17) base k
        vDF = vDF*score_RNA(10)*score_RNA(29+nphitype(JN))

!     END ENERGY WITH RESPECT TO COSPHI


!     ----- DC = FIRST DER. OF COSPHI W/RESPECT TO THE CARTESIAN DIFFERENCES T -----
        rDC = -rG*Z12-vCPHI*rD*Z10**2
        rDC2 = rD*Z12+vCPHI*rG*Z20**2
!     ----- UPDATE THE FIRST DERIVATIVE ARRAY -----
        rDR1 = vDF*(crossproduct(rKJ,rDC))
        rDR2 = vDF*(crossproduct(rKJ,rDC2))
        rDR = vDF*(crossproduct(rIJ,rDC) + crossproduct(rDC2, rKL))
        rFI = - rDR1
        rFJ = - rDR + rDR1
        rFK = + rDR + rDR2
        rFL = - rDR2

        F(I3+1:I3+3) = F(I3+1:I3+3) + rFI
        F(J3+1:J3+3) = F(J3+1:J3+3) + rFJ
        F(K3+1:K3+3) = F(K3+1:K3+3) + rFK
        F(L3+1:L3+3) = F(L3+1:L3+3) + rFL
        ETORS = ETORS + vEPW
     END DO
     
  END SUBROUTINE ENERGY_DIHS

   SUBROUTINE DEALLOC_DIHS()
      IF (ALLOCATED(PK)) DEALLOCATE(PK)
      IF (ALLOCATED(PN)) DEALLOCATE(PN)
      IF (ALLOCATED(IPN)) DEALLOCATE(IPN)
      IF (ALLOCATED(PHASE)) DEALLOCATE(PHASE)
      IF (ALLOCATED(IP)) DEALLOCATE(IP)
      IF (ALLOCATED(JP)) DEALLOCATE(JP)
      IF (ALLOCATED(KP)) DEALLOCATE(KP)
      IF (ALLOCATED(LP)) DEALLOCATE(LP)
      IF (ALLOCATED(ICP)) DEALLOCATE(ICP)
      IF (ALLOCATED(GAMC)) DEALLOCATE(GAMC)
      IF (ALLOCATED(GAMS)) DEALLOCATE(GAMS)
      IF (ALLOCATED(FMN)) DEALLOCATE(FMN)
      IF (ALLOCATED(NPHITYPE)) DEALLOCATE(NPHITYPE)
   END SUBROUTINE DEALLOC_DIHS

END MODULE MOD_DIHEDRALS






