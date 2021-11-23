MODULE MOD_EXCLV
  USE PREC_HIRE
  USE NBDEFS
  IMPLICIT NONE
  REAL(KIND = REAL64) :: EXCL_VOL
  REAL(KIND = REAL64) :: BARRIER
  REAL(KIND = REAL64) :: RATIO
  
  CONTAINS

  ! Setting the common variables to the correct values
  SUBROUTINE INIT_EXCLV()
     USE NAPARAMS, ONLY: SCORE_RNA
     
     !TODO: When extending to proteins make these arrays for each residue?
     EXCL_VOL = SCORE_RNA(11)  ! steepness
     BARRIER = SCORE_RNA(44)  ! Height = 100
     RATIO = SCORE_RNA(45)    ! Not used with old barrier
  
  END SUBROUTINE INIT_EXCLV

  SUBROUTINE ENERGY_EXCLV(DA2, EEXCL, DF, CT2, R2IN, R2OUT)   

     REAL(KIND = REAL64), INTENT(IN) :: DA2        !input squared distance between particles
     REAL(KIND = REAL64), INTENT(IN) :: CT2        !QUERY: what is this?
     REAL(KIND = REAL64), INTENT(IN) :: R2IN       !lower cutoff
     REAL(KIND = REAL64), INTENT(IN) :: R2OUT      !higher cutoff
     REAL(KIND = REAL64), INTENT(OUT) :: DF        !force contribution
     REAL(KIND = REAL64), INTENT(OUT) :: EEXCL     !energy contribution

     REAL(KIND = REAL64) :: R, CT2MOD, EXPEXCL


     R = DSQRT(DA2)
     CT2MOD = RATIO*CT2 ! introduced by sp Apr20
     ! CT2MOD = 1.0D0/BARRIER*(LOG(0.8D0/(BARRIER-0.380)))+CT2 ! old barrier
     EXPEXCL = EXP(-EXCL_VOL*(R-CT2MOD))
     EEXCL = BARRIER*(1-1/(1+EXPEXCL))
     DF = -BARRIER*EXCL_VOL*EXPEXCL/(R*(1+EXPEXCL)**2)       
!     Eexcl = excl_vol*exp(4.0*(ct2-r))
!     DF = -excl_vol*(4.0*exp(4.0*(ct2-r)))/r 
     IF (DA2 .GE. R2IN) CALL RNA_SWITCH_CUTOFF(DA2,EEXCL,DF,R2IN,R2OUT)
  END SUBROUTINE ENERGY_EXCLV
  
  SUBROUTINE RNA_SWITCH_CUTOFF(R2, ESW, FSW, RI2, RO2)
    REAL(KIND = REAL64), INTENT(IN) :: R2, RI2, RO2
    REAL(KIND = REAL64), INTENT(INOUT) :: ESW  !energy
    REAL(KIND = REAL64), INTENT(INOUT) :: FSW  !force
    
    REAL(KIND = REAL64) :: RD6, SW, DSW
    
    RD6 = 1.0D0/(RO2-RI2)**3
    !QUERY: Should we set these magic numbers as parameters?
    SW = (RO2+2.0D0*R2-3.0D0*RI2)*RD6*(RO2-R2)**2
    DSW = 12.0D0*(RO2-R2)*(RI2-R2)*RD6
    
    FSW = FSW*SW - ESW*DSW
    ESW = ESW*SW
  END SUBROUTINE RNA_SWITCH_CUTOFF
  
END MODULE MOD_EXCLV



