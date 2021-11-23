MODULE MOD_BASESTACKING
  USE PREC_HIRE
  USE NBDEFS
  IMPLICIT NONE

  
  CONTAINS
  
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

! Stacking interactions
!       i-2    i   j-2    j
!         \   /      \   /
!          \ /        \ /
!          i-1        j-1
  
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
END MODULE MOD_BASESTACKING
