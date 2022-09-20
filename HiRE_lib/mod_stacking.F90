!> @file
!> Contains MOD_BASESTACKING module to provide base stacking contributions

!> Setup and calculation of base stacking contributions.
!> Two version exist - the original implementation and one with a vertical offset.
MODULE MOD_BASESTACKING
   USE PREC_HIRE
   USE NBDEFS
   IMPLICIT NONE

   CONTAINS
      !> Old stacking parameters
      !> @brief
      !> 
      !> Routine returns the stacking parameters for a given base pair.\n
      !> The subroutine differentiates between pyr-pyr, pur-pur and pyr-pur interactions.
      !>
      !> @param[in] TI - type of base I
      !> @param[in] TJ - type of base J
      !> @param[out] EQ - stacking equilibrium distance
      !> @param[out] WID - stacking width parameter (lateral variance that still allows interactions)   
      !> @param[out] SK - scaling of stacking interactions
      SUBROUTINE STACKPARAMS(TI,TJ,EQ,WID,SK)
         USE NAPARAMS, ONLY: SCORE_RNA
         INTEGER, INTENT(IN) :: TI, TJ
         REAL(KIND = REAL64), INTENT(OUT) :: SK, EQ, WID
            
         IF ((TI.LT.3.AND.TJ.GT.2).OR.(TJ.LT.3.AND.TI.GT.2)) THEN
            !pyr-pur
            EQ = SCORE_RNA(17)
            WID = SCORE_RNA(20)
            SK = SCORE_RNA(14)
         ELSE IF (TI.LT.3.AND.TJ.LT.3) THEN
            !pur-pur
            EQ = SCORE_RNA(18)
            WID = SCORE_RNA(21)
            SK = SCORE_RNA(15)
         ELSE
            !pyr-pyr
            EQ = SCORE_RNA(19)
            WID = SCORE_RNA(22)
            SK = SCORE_RNA(16)       
         ENDIF           
      END SUBROUTINE STACKPARAMS

      !> New stacking parameters
      !> @brief
      !> 
      !> Routine returns the stacking parameters for a given base pair.\n
      !> The subroutine differentiates between pyr-pyr, pur-pur and pyr-pur interactions.
      !>
      !> @param[in] TI - type of base I
      !> @param[in] TJ - type of base J
      !> @param[out] EQ - stacking equilibrium distance
      !> @param[out] WID - stacking width parameter (lateral variance that still allows interactions)   
      !> @param[out] SK - scaling of stacking interactions
      !> @param[out] Th - angular part for offset
      !> @param[out] GM - Query - what is this?          
      SUBROUTINE STACKPARAMS2(TI,TJ,EQ,WID,SK,Th,GM)
         USE NAPARAMS, ONLY: SCORE_RNA
         INTEGER, INTENT(IN) :: TI, TJ
         REAL(KIND = REAL64), INTENT(OUT) :: SK, EQ, WID, Th, GM
            
         IF ((TI.LT.3.AND.TJ.GT.2).OR.(TJ.LT.3.AND.TI.GT.2)) THEN
            !pyr-pur
            EQ = SCORE_RNA(17)
            WID = SCORE_RNA(20)
            SK = SCORE_RNA(14)
            Th = 20 * 3.15159/180             ! hard coded for testing --> to be added to SCORE_RNA
            GM = 8                            ! hard coded for testing --> to be added to SCORE_RNA
      !          Th = SCORE(48)
      !          GM = SCORE(49)
         ELSE IF (TI.LT.3.AND.TJ.LT.3) THEN
            !pur-pur
            EQ = SCORE_RNA(18)
            WID = SCORE_RNA(21)
            SK = SCORE_RNA(15)
            Th = 35 * 3.15159/180          ! hard coded for testing --> to be added to SCORE_RNA
            GM = 8                            ! hard coded for testing --> to be added to SCORE_RNA
      !          Th = SCORE(50)
      !          GM = SCORE(51)
         ELSE
            !pyr-pyr
            EQ = SCORE_RNA(19)
            WID = SCORE_RNA(22)
            SK = SCORE_RNA(16)   
            Th = 40 * 3.15159/180             ! hard coded for testing --> to be added to SCORE_RNA
            GM = 8                            ! hard coded for testing --> to be added to SCORE_RNA
      !          Th = SCORE(52)
      !          GM = SCORE(53)
         ENDIF       
      END SUBROUTINE STACKPARAMS2
    

    
! Stacking interactions
!       i-2    i   j-2    j
!         \   /      \   /
!          \ /        \ /
!          i-1        j-1
 
      !> Old stacking potential
      !> @brief
      !>
      !> Old energy and gradient for the base stacking potential, without a vertical offset.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] I - index of final particle in base 1
      !> @param[in] J - index of final particle in base 2
      !> @param[in] TI - type of base 1
      !> @param[in] TJ - type of base 2
      !> @param[in] X - input coordinates
      !> @param[out] ESTK - stacking energy for the given base pair
      !> @param[out] F - gradient from stacking interactions
      !> @param[in] STACKUNIT - output unit for debugging
    SUBROUTINE RNA_STACKV(NOPT,I,J,TI,TJ,F,X,ESTK,STACKUNIT)
      USE VEC_UTILS
      USE NAPARAMS, ONLY: BTYPE
      INTEGER, INTENT(IN) :: NOPT                   !Number of degrees of freedom
      INTEGER, INTENT(IN) :: I, J                   !indices of final particle in residue
      INTEGER, INTENT(IN) :: TI,TJ
      INTEGER, INTENT(IN) :: STACKUNIT              !output unit, only assigned to file if FOR_ANALYSIS is used as compile flag
      REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
      REAL(KIND = REAL64), INTENT(INOUT) :: F(NOPT)   !force
      REAL(KIND = REAL64), INTENT(OUT) :: ESTK

      REAL(KIND = REAL64), PARAMETER :: COSP=2.0 !Power of cos(x) in energy function
      REAL(KIND = REAL64) :: A(3), B(3), C(3), D(3)  !Particle coordinates
      REAL(KIND = REAL64) :: AxB0(3), CxD0(3), VA, VC
      REAL(KIND = REAL64) :: R(3), R0(3), R1, R2, DotP, dot1, dot2
      REAL(KIND = REAL64) :: dot1w, dot2w, dot1wd, dot2wd, dotw, dotwd
      REAL(KIND = REAL64) :: Dvr(3), Dvrij(3), Da(3), Fx, Fy, Fz
      REAL(KIND = REAL64) :: SK, EQ, WID     
      REAL(KIND = REAL64), PARAMETER :: EPS = 1.0D-6

      CALL STACKPARAMS(TI,TJ,EQ,WID,SK)
      
      !Get relevant coordinates
      A(1:3) = X((3*(I-2)-2):(3*(I-2))) - X((3*(I-1)-2):(3*(I-1))) !! vector a : I-1 -> I-2
      B(1:3) = X((3*I-2):(3*I))         - X((3*(I-1)-2):(3*(I-1))) !! vector b : I-1 -> I
      C(1:3) = X((3*(J-2)-2):(3*(J-2))) - X((3*(J-1)-2):(3*(J-1))) !! vector c : J-1 -> J-2
      D(1:3) = X((3*J-2):(3*J))         - X((3*(J-1)-2):(3*(J-1))) !! vector d : J-1 -> J
      
      !Get normalised crossproducts
      CALL NORMED_CP(A,B, AxB0, VA)
      CALL NORMED_CP(C,D, CxD0, VC)
      
      !Get vector R
      R(1:3)= (X((3*(I-2)-2):(3*(I-2))) + X((3*(I-1)-2):(3*(I-1))) &
             + X((3*I-2):(3*I))         - X((3*(J-2)-2):(3*(J-2)))  &
             - X((3*(J-1)-2):(3*(J-1))) - X((3*J-2):(3*J)))/3
      R1 = EUC_NORM(R)
      R0(1:3) = (1.0D0/R1)*R(1:3)   
      
      !Various dot products of different vectors
      DotP = dot_product(axb0, cxd0)
      dotw = 1 - (1 -DotP**2)**2                ! |ni x nj|^4
      dotwd = 2*(1-DotP**2) / (2-dotP**2)
      dot1 = dot_product(r0, axb0)              
      dot1w = 1 - (1 -dot1**2)**2               ! |ni x r|^4
      dot1wd = 2*(1-dot1**2) / (2-dot1**2)
      dot2 = dot_product(r0, cxd0)             
      dot2w = 1 - (1 -dot2**2)**2               ! |nj x r|^4
      dot2wd = 2*(1-dot2**2) / (2-dot2**2)      
      
      R2 = (R1 - EQ)/WID
      !QUERY: So which one is it?
!      Estk = -SK * DotP**cosp * dexp(-r2**2)
!      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1**cosp * dot2**cosp
!      Estk = -SK * DotP**cosp * dexp(-r2**2) * dot1w * dot2w
      Estk = -SK * dotw * dexp(-r2**2) * dot1w * dot2w           !  ciccata ????       

#if FOR_ANALYSIS
      IF (ABS(ESTK) .GT. 0.1D0) THEN
            WRITE(STACKUNIT,'(4I6,14F15.7)') I, J, ti, tj, -SK, r1, eq, wid, dexp(-r2**2), DotP,DotP**cosP, dot1, &
            (1 -dot1**2)**2, dot1w, dot2, (1 -dot2**2)**2, dot1w, Estk
      ENDIF
#endif

      ! Bypass derivatives calculation if E is very small
      ! This also prevents unstabilities arising from cos(x) ~= 0      
      IF ((ESTK.GT.(-EPS)) .OR. (ABS(DOT1).LT.EPS) .OR. (ABS(DOT2).LT.EPS) &
          .OR. (ABS(DOTP).LT.EPS)) THEN
          ESTK = 0.0D0
          RETURN
      ENDIF
    
!      Estk = -SK * DotP**4 * (eq/r1)**6
!      Estk = -SK * DotP**4 * dexp(-3*(r1-3))
      Dvr = -Estk * 1/3 * 2*r2**1/wid * r/r1
!      Dvr = -Estk * 1/2 * 6/r1  *r/r1
!      Dvr = -Estk * 1/3 * (-3*(r1-3))  * -3*r/r1
!      Dvrij = Estk*cosp * 1/(3*r1) * (axb0/dot1 + cxd0/dot2 - 2*r0)
      Dvrij = Estk*cosp * 1/(3*r1) * (dot1wd*axb0/dot1 + dot2wd*cxd0/dot2 - &
              r0*(dot1wd+dot2wd))
      Dvr = Dvr + Dvrij

!------- Derivatives on the 6 particles   --

!      Da = (cxd0/DotP - axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + r0/dot1 - 2*axb0)*cosp*Estk/VA
!      Da = (cxd0/DotP + dot1wd*r0/dot1 - axb0*(1+dot1wd))*cosp*Estk/VA
      Da = (dotwd*cxd0/DotP + dot1wd*r0/dot1 - axb0*(dotwd+dot1wd))*2*Estk/VA
!      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - crossproduct(b, Da)
      F(i*3-8:i*3-6) = F(i*3-8:i*3-6) - Dvr - crossproduct(b, Da)
      F(i*3-5:i*3-3) = F(i*3-5:i*3-3) - Dvr - crossproduct(a-b, Da)
      F(i*3-2:i*3  ) = F(i*3-2:i*3  ) - Dvr - crossproduct(Da, a)
!      Da = (axb0/DotP - cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + r0/dot2 - 2*cxd0)*cosp*Estk/VC
!      Da = (axb0/DotP + dot2wd*r0/dot2 - cxd0*(1+dot2wd))*cosp*Estk/VC
      Da = (dotwd*axb0/DotP + dot2wd*r0/dot2 - cxd0*(dotwd+dot2wd))*2*Estk/VC
!      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) - crossproduct(d, Da)
      F(j*3-8:j*3-6) = F(j*3-8:j*3-6) + Dvr - crossproduct(d, Da)
      F(j*3-5:j*3-3) = F(j*3-5:j*3-3) + Dvr - crossproduct(c-d, Da)
      F(j*3-2:j*3  ) = F(j*3-2:j*3  ) + Dvr - crossproduct(Da, c)
    END SUBROUTINE RNA_STACKV 
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! NEW STAKING POTENTIAL - Vertical offset !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
  
      !> New stacking potential
      !> @brief
      !>
      !> New energy and gradient for the base stacking potential, including a vertical offset.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] I - index of final particle in base 1
      !> @param[in] J - index of final particle in base 2
      !> @param[in] TI - type of base 1
      !> @param[in] TJ - type of base 2
      !> @param[in] X - input coordinates
      !> @param[out] ESTK - stacking energy for the given base pair
      !> @param[out] F - gradient from stacking interactions
      SUBROUTINE RNA_STACKV2(NOPT,I,J,TI,TJ,F,X,ESTK)
         USE VEC_UTILS
         USE NAPARAMS, ONLY: BTYPE
         INTEGER, INTENT(IN) :: NOPT                   !Number of degrees of freedom
         INTEGER, INTENT(IN) :: I, J                   !indices of final particle in residue
         INTEGER, INTENT(IN) :: TI,TJ
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force
         REAL(KIND = REAL64), INTENT(OUT) :: ESTK
      
         REAL(KIND = REAL64), PARAMETER :: COSP=2.0 !Power of cos(x) in energy function
         REAL(KIND = REAL64) :: A(3), B(3), C(3), D(3)  !Particle coordinates
         REAL(KIND = REAL64) :: AxB0(3), CxD0(3), AxB(3), CxD(3), VA, VC
         REAL(KIND = REAL64) :: R(3), R0(3), R1, R2, DotP
         REAL(KIND = REAL64) :: ct01, ct02, st01, st02, CosT, CosT01, CosT02
         REAL(KIND = REAL64) :: Estkrs, Dva1, Dva2, Dvrg(3), DvraT(3), Dsp
         REAL(KIND = REAL64) :: adb, adc, add, bdc, bdd, cdd, a2, b2, c2, d2, dxabcd, drcxd
         REAL(KIND = REAL64) :: rxd(3), rxc(3)
         REAL(KIND = REAL64) :: DIs(3), DIm1s(3), DIm2s(3), DJs(3), DJm1s(3), DJm2s(3), DJv(3), DJm1v(3), DJm2v(3)
         REAL(KIND = REAL64) :: Dvr(3), Dvrij(3), Da(3), Fx, Fy, Fz
         REAL(KIND = REAL64) :: SK, EQ, WID, Th, GM     
         REAL(KIND = REAL64), PARAMETER :: EPS = 1.0D-6

         CALL STACKPARAMS2(TI,TJ,EQ,WID,SK,Th,GM)          ! two new parameters Th and GM to add 
      
         !Get relevant coordinates
         A(1:3) = X((3*(I-2)-2):(3*(I-2))) - X((3*(I-1)-2):(3*(I-1))) !! vector a : I-1 -> I-2
         B(1:3) = X((3*I-2):(3*I))         - X((3*(I-1)-2):(3*(I-1))) !! vector b : I-1 -> I
         C(1:3) = X((3*(J-2)-2):(3*(J-2))) - X((3*(J-1)-2):(3*(J-1))) !! vector c : J-1 -> J-2
         D(1:3) = X((3*J-2):(3*J))         - X((3*(J-1)-2):(3*(J-1))) !! vector d : J-1 -> J
      
         !Get normalised crossproducts   --> give spin vectors
         CALL NORMED_CP2(A,B, AXB, AxB0, VA)
         CALL NORMED_CP2(C,D, CXD, CxD0, VC)
         DotP = dot_product(axb0, cxd0)                ! spins orientation

         !Get vector R  --> vector connecting the bases
         R(1:3)= (X((3*(I-2)-2):(3*(I-2))) + X((3*(I-1)-2):(3*(I-1))) &
                + X((3*I-2):(3*I))         - X((3*(J-2)-2):(3*(J-2)))  &   
                - X((3*(J-1)-2):(3*(J-1))) - X((3*J-2):(3*J)))/3
         R1 = EUC_NORM(R)
         R0(1:3) = (1.0D0/R1)*R(1:3)       
         R2 = (R1 - EQ)/WID                             ! base distance computed from the center of mass of the 3 particles
      
         !Get vertical angular position cos(theta - theta_0) --> CosTT0
!         Th=0
         ct01 = cos(Th)
         st01 = sin(Th)
         ct02 = cos(3.14-Th)                            ! Pi - Theta_0
         st02 = sin(3.14-Th)
         CosT = dot_product(R0, cxd0)
         CosT01 = CosT*ct01 + sqrt(1-CosT*CosT)*st01
         CosT02 = CosT*ct02 + sqrt(1-CosT*CosT)*st02

         Estkrs = -SK * dexp(-r2**2) * DotP * DotP
         Estk = Estkrs*(exp(-GM*(1-CosT01))+exp(-GM*(1-CosT02)))
!         Estkrs = 1
      
         Dsp = -2*SK * dexp(-r2**2)  * DotP *(exp(-GM*(1-CosT01))+exp(-GM*(1-CosT02)))

      ! Bypass derivatives calculation if E is very small
      ! This also prevents unstabilities arising from cos(x) ~= 0      
!       IF ((ESTK.GT.(-EPS)) .OR. (ABS(DOTP).LT.EPS) .OR. (dexp(-r2**2).LT.EPS)) THEN
!           ESTK = 0.0D0
!           RETURN
!       ENDIF
      
!       IF (ABS(ESTK) .GT. 0.0D0) THEN
!             WRITE(*,'(4I6,14F15.7)') I, J, ti, tj,  CosT01, CosT02, dexp(-GM*(1-CosT01)), dexp(-GM*(1-CosT02)), Estk
!       ENDIF

         ! Global derivatives
         Dva1 = GM*(ct01 - st01*CosT/sqrt(1-CosT*CosT))*exp(-GM*(1-CosT01))   ! derivative of vertical exponential offset  
         Dva2 = GM*(ct02 - st02*CosT/sqrt(1-CosT*CosT))*exp(-GM*(1-CosT02))   ! derivative of vertical exponential offset 
      
         Dvrg = -Estk * 2*r2**1/wid * R/R1                                                 ! derivative of gaussian over distance over R
         DvraT = (Dva1 + Dva2)*Estkrs*(cxd0 - CosT*R0)/R1                                  ! derivative of vertical exponential offset over R
         Dvr = (Dvrg + DvraT)/3                                                            ! total derivative over R

         ! dot and cross products for derivatives
         adb = dot_product(A,B)
         adc = dot_product(A,C)
         add = dot_product(A,D)
         bdc = dot_product(B,C)
         bdd = dot_product(B,D)
         cdd = dot_product(C,D)
         a2 = dot_product(A,A)
         b2 = dot_product(B,B)
         c2 = dot_product(C,C)
         d2 = dot_product(D,D)
         dxabcd = dot_product(axb,cxd)                      ! not necessary, use DotP instead
         rxc = crossproduct(R0,C)/VC
         rxd = crossproduct(R0,D)/VC

         ! Derivatives on the 6 particles of the spins contribution
         DIs = (adc*D - add*C - dxabcd*(a2*B - adb*A)/(VA*VA))/(VA*VC)
         DIm1s = (dot_product(D,A-B)*C + dot_product(C,B-A)*D - dxabcd*(adb*(A+B) - a2*B - b2*A)/(VA*VA))/(VA*VC)
         DIm2s = (bdd*C - bdc*D- dxabcd*(b2*A- adb*B)/(VA*VA))/(VA*VC)
         DJs = (adc*B - bdc*A - dxabcd*(c2*D - cdd*C)/(VC*VC))/(VA*VC)
         DJm1s = (dot_product(B,C-D)*A + dot_product(A,D-C)*B - dxabcd*(cdd*(C+D) - c2*D - d2*C)/(VC*VC))/(VA*VC)
         DJm2s = (bdd*A - add*B- dxabcd*(d2*C- cdd*D)/(VC*VC))/(VA*VC)
      
         ! Derivatives on the 3 J particles of the vertical contribution
         DJv = rxc - CosT*(c2*D - cdd*C)/(VC*VC)
         DJm1v = rxd - rxc - CosT*(cdd*(C+D) - d2*C - c2*D)/(VC*VC)
         DJm2v = -rxd - CosT*(d2*C - cdd*D)/(VC*VC)
      
         F(i*3-8:i*3-6) = F(i*3-8:i*3-6) -Dvr - Dsp*DIm2s 
         F(i*3-5:i*3-3) = F(i*3-5:i*3-3) -Dvr - Dsp*DIm1s
         F(i*3-2:i*3  ) = F(i*3-2:i*3  ) -Dvr - Dsp*DIs

         F(j*3-8:j*3-6) = F(j*3-8:j*3-6) + Dvr - Dsp*DJm2s - Estkrs*(Dva1+Dva2)*DJm2v
         F(j*3-5:j*3-3) = F(j*3-5:j*3-3) + Dvr - Dsp*DJm1s - Estkrs*(Dva1+Dva2)*DJm1v
         F(j*3-2:j*3  ) = F(j*3-2:j*3  ) + Dvr - Dsp*DJs   - Estkrs*(Dva1+Dva2)*DJv
      
      
      END SUBROUTINE RNA_STACKV2 
    
END MODULE MOD_BASESTACKING
