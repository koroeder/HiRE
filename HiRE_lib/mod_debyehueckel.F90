MODULE MOD_DEBYEHUECKEL
  USE PREC_HIRE
  USE NBDEFS
  IMPLICIT NONE
  REAL(KIND = REAL64) :: DL  
  REAL(KIND = REAL64) :: DIEL 
    
  CONTAINS

  SUBROUTINE INIT_DH()
     USE NAPARAMS, ONLY: SCORE_RNA
     DIEL = SCORE_RNA(12)
     DL = SCORE_RNA(13)  
  END SUBROUTINE INIT_DH

  SUBROUTINE ENERGY_DH(NOPT, X, F, EDH)
     USE VAR_DEFS, ONLY: CHATM, NRES, RESSTART, RESFINAL
     USE UTILS_IO, ONLY: GETUNIT
     IMPLICIT NONE
     
     INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
     REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
     REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
     REAL(KIND = REAL64), INTENT(OUT) :: EDH

     INTEGER :: I, J, K, L
     REAL(KIND = REAL64) :: RIJ(3), R, CHRGI, CHRGJ, EDHPAIR, DFPAIR, NB, DX(3)

#if FOR_ANALYSIS
     INTEGER :: DHUNIT

     DHUNIT = GETUNIT()
     OPEN(DHUNIT, FILE="Dbg_DebyeHueckel.dat", STATUS='UNKNOWN')
#endif

     EDH = 0.0D0
     F(1:NOPT) = 0.0D0
     
     !QUERY: Do we want all charge-charge interactions, or only does in between
     !       atoms in different residues?
     !TODO: this loop needs to be over the residues and then we iterate inside over the particles!
     DO K=1,NRES-1
       DO L=K+1,NRES
         DO I = RESSTART(K),RESFINAL(K)
           DO J = RESSTART(L),RESFINAL(L)
             RIJ(1:3) = X(I*3-2:I*3) - X(3*J-2:3*J)
             R = DSQRT(DOT_PRODUCT(RIJ,RIJ))
             CHRGI = CHATM(I)
             CHRGJ = CHATM(J)
             EDHPAIR = 0.0D0
             DFPAIR = 0.0D0
             ! if charges are non-zero, calculate Debye-Hueckel contribution
             IF (ABS(CHRGI*CHRGJ).GT.1.0D-6) THEN
                CALL DH_PAIR(R, EDHPAIR, DFPAIR, CHRGI, CHRGJ)
                NB = GET_NBCOEF(I,J)
                EDH = EDH + EDHPAIR * NB               
                DX(1:3) = DFPAIR * NB * RIJ(1:3)
#ifdef FOR_ANALYSIS
                WRITE(DHUNIT,'(4I6,2F7.3,4F15.7)') K, L, I, J, CHRGI, CHRGJ, R, NB, EDHPAIR, EDHPAIR*NB
#endif
                F((I*3-2):I*3) = F((I*3-2):I*3) - DX(1:3)
                F((J*3-2):J*3) = F((J*3-2):J*3) + DX(1:3) 
             ENDIF
           END DO
         END DO
       END DO
     END DO         
#ifdef FOR_ANALYSIS
     CLOSE(DHUNIT)
#endif
  END SUBROUTINE ENERGY_DH
  
  SUBROUTINE DH_PAIR(R, EDH, DF, QI, QJ)
     REAL(KIND = REAL64), INTENT(IN) :: R      ! distance between i and j
     REAL(KIND = REAL64), INTENT(IN) :: QI, QJ ! charges on i and j
     REAL(KIND = REAL64), INTENT(OUT) :: EDH   ! energy for this pair
     REAL(KIND = REAL64), INTENT(OUT) :: DF    ! force for this pair   

     ! Debye-Huckel prefactor : 1/(4*pi*e_0*e_r)
     !   1/(4*pi*e_0) ~= 332.056 kcal/mol * A
     !   1/(4*pi*e_0*e_r) ~= 332.056/80 = 4.1507
     ! Debye length = 1/sqrt(8*pi*l_b*I)
     !   l_b ~= 7
     !   I : Ionic strength
     !   DL ~= 3.04 / sqrt(I)    
     EDH = DIEL*(4.1507*QI*QJ/R)*EXP(-R/DL)
     DF = -EDH*(1/DL+1/R)/R
     !write(83,*) qi, qj, r, "  ", diel, dl, Edh, DF 
  END SUBROUTINE DH_PAIR
  
END MODULE MOD_DEBYEHUECKEL

