!> @file
!> Contains MOD_DIHEDRALS to handle dihedral angles

!> Module containing all routines and variables to calculate the dihedral contribution to E and F\n
!> Dihedrals: E = PK * (1 + COS(PN * P - PHASE))
MODULE MOD_DIHEDRALS
   USE PREC_HIRE
   !> Number of dihedrals
   INTEGER :: NDIHS   
   !> Number of unique dihedral types
   INTEGER :: NPTRA    !Number of unique dihedral types

   !Dihedrals: E = PK * (1 + COS(PN * P - PHASE)) with P the current dihedral angle
   !> Torsional force constant
   REAL(KIND = REAL64), ALLOCATABLE :: PK(:,:)    
   !> Number of torsional equilibrium states 
   REAL(KIND = REAL64), ALLOCATABLE :: PN(:,:)
   !> Interger value of PN    
   INTEGER, ALLOCATABLE :: IPN(:)               
   !> Torsional phase factor
   REAL(KIND = REAL64), ALLOCATABLE :: PHASE(:,:)
   !> Torsional y-offset
   REAL(KIND = REAL64), ALLOCATABLE :: DOFFSET(:)  
   !> Reference values for torsions (cosine)
   REAL(KIND = REAL64), ALLOCATABLE :: GAMC(:,:)  
   !> Reference values for torsions (sine)
   REAL(KIND = REAL64), ALLOCATABLE :: GAMS(:,:)   

   !book-keeping for torsionals
   !> Dihedral atom1
   INTEGER, ALLOCATABLE :: IP(:)                 
   !> Dihedral atom2
   INTEGER, ALLOCATABLE :: JP(:)                 
   !> Dihedral atom3
   INTEGER, ALLOCATABLE :: KP(:)                 
   !> Dihedral atom4
   INTEGER, ALLOCATABLE :: LP(:)                 
   !> Type of dihedral
   INTEGER, ALLOCATABLE :: ICP(:)                

   CONTAINS
      !> Routine to allocate all required arrays
      SUBROUTINE ALLOC_DIHS()
         CALL DEALLOC_DIHS()
         ALLOCATE(PK(NPTRA,3), PN(NPTRA,3), IPN(NPTRA), PHASE(NPTRA,3), &
                  IP(NDIHS), JP(NDIHS), KP(NDIHS), LP(NDIHS), ICP(NDIHS), &
                  GAMC(NPTRA,3), GAMS(NPTRA,3), DOFFSET(NPTRA))
      END SUBROUTINE ALLOC_DIHS 

      !> Initialise helper variables to deal with phase shifting etc.
      SUBROUTINE INIT_DIHPAR()
         USE NUM_DEFS, ONLY: PI
         IMPLICIT NONE
         REAL(KIND = REAL64) :: DUM, DUMS, DUMC
         REAL(KIND = REAL64), PARAMETER :: EPS1 = 1.0D-3
         REAL(KIND = REAL64), PARAMETER :: EPS2 = 1.0D-6
         INTEGER :: I,J

         DO I = 1,NPTRA
            DO J=1,3
               IF (PK(I,J).GT.0.0D0) THEN
                  DUM = PHASE(I,J)
                  IF (DABS(DUM-PI) .LE. EPS1) DUM = SIGN(PI,DUM)
                  DUMC = DCOS(DUM)
                  DUMS = DSIN(DUM)
                  IF(DABS(DUMC) .LE. EPS2) DUMC = 0.0d0
                  IF(DABS(DUMS) .LE. EPS2) DUMS = 0.0d0

                  GAMC(I,J) = DUMC*PK(I,J)
                  GAMS(I,J) = DUMS*PK(I,J)
         
                  PN(I,J) = DABS(PN(I,J))
                  IPN(I) = INT(PN(I,J)+EPS1)  
                  IF (IPN(I).EQ.0) IPN(I) = 1 !needed for zero terms, GMUl of 1 is 0.0, sothis does not change the energy computation but prevents a seg fault
               END IF
            END DO    
         ENDDO
      END SUBROUTINE INIT_DIHPAR

      !> Routine to pass dihedral information to external programmes
      SUBROUTINE GET_DIHEDRALS(DIHINFO)
         INTEGER, INTENT(OUT) :: DIHINFO(NDIHS,4)
         INTEGER JN

         DO JN=1,NDIHS
            DIHINFO(JN,1) = IP(JN)/3 + 1
            DIHINFO(JN,2) = JP(JN)/3 + 1
            DIHINFO(JN,3) = IABS(KP(JN))/3 + 1
            DIHINFO(JN,4) = IABS(LP(JN))/3 + 1      
         ENDDO
      END SUBROUTINE GET_DIHEDRALS

      !> Calculate the energy and force contribution from the dihedral terms
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
         REAL(KIND = REAL64), PARAMETER :: EPS0 = 1.0D-5
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

         INTEGER :: JN, IC, INC, K
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
            !iterate over all terms
            DO K=1,3
               IC = ICP(JN)
               INC = IPN(IC)
               IF (ABS(PK(IC,K)).GT.EPS0) THEN
                  CT0 = PN(IC,K)*AP1
                  COSNP = DCOS(CT0)
                  SINNP = DSIN(CT0)

                  !not sure the below is needed?
                  !if (PN(IC).eq.12) then
                  !   !GAMCs=GAMC(IC)/PK(IC)
                  !   !GAMSs=GAMS(IC)/PK(IC)
                  !   COSNPs=DCOS(AP1)
                  !   SINNPs=DSIN(AP1)
                  !   EXPRB=ACOS(GAMC(IC)/PK(IC))*180.d0/PI
         
                  !   vEPW=(PK(IC)*COSNPs**int(EXPRB))*vFMUL

                  !   if (EXPRB.eq.0) then
                  !      DF0=0.d0
                  !      df1=0.d0
                  !   else
                  !      DF0=PK(IC)*EXPRB*SINNPs*COSNPs**int(EXPRB-1.d0)
                  !      !print*,DF0
                  !      DUMS = vSPHI+SIGN(TM24,vSPHI)
                  !      DFLIM = GAMC(IC)*(PN(IC)-GMUL(INC)+GMUL(INC)*vCPHI)
                  !   
                  !      df1 = df0/dums
                  !      if (tm06.gt.abs(dums)) df1 = dflim
                  !   endif
                  !   !DF0 = -PN(IC)*PK(IC)*(GAMCs*SINNP-GAMSs*COSNP)*(COSNP*GAMCs+SINNP*GAMSs)**(PN(IC)-1.d0)
               !else
                     vEPW= (PK(IC,K)+COSNP*GAMC(IC,K)+SINNP*GAMS(IC,K))*vFMUL !! might be revised
                     DF0 = PN(IC,K)*(GAMC(IC,K)*SINNP-GAMS(IC,K)*COSNP)
                     DUMS = vSPHI+SIGN(TM24,vSPHI)
                     DFLIM = GAMC(IC,K)*(PN(IC,K)-GMUL(INC)+GMUL(INC)*vCPHI)
                        
                     df1 = df0/dums
                     if(tm06.gt.abs(dums)) df1 = dflim
                  !endif
                  
                  vDF = DF1*vFMUL

                  vEPW = vEPW*SCORE_RNA(4) 
                  vDF = vDF*SCORE_RNA(4)

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
               END IF
            END DO
         END DO   
      END SUBROUTINE ENERGY_DIHS

      !> Deallocate all arrays in this module
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
      END SUBROUTINE DEALLOC_DIHS

END MODULE MOD_DIHEDRALS






