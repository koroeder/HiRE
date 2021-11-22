MODULE CALCFORCES
  USE PREC_HIRE
  IMPLICIT NONE

  TYPE POT_ENE
     REAL(KIND = REAL64) :: EBOND         !bonding energy
     REAL(KIND = REAL64) :: EANGLES       !bond angles
     REAL(KIND = REAL64) :: ETORS         !torsional energy
     REAL(KIND = REAL64) :: EDH           !Debye-Hueckel energy     
     REAL(KIND = REAL64) :: EHBOND        !hydrogen bonding
     REAL(KIND = REAL64) :: EVDW          !excluded volume 
     REAL(KIND = REAL64) :: ESTAK         !stacking interactions
     REAL(KIND = REAL64) :: ESAXS         !SAXS energy 
     REAL(KIND = REAL64) :: EDISTR        !Distance contraints
     REAL(KIND = REAL64) :: EPOSR         !Positional contraints 
     REAL(KIND = REAL64) :: ETOT          !Total energy 
  END TYPE POT_ENE
  
  TYPE(POT_ENE) :: EVEC  !vector containing the energy contributions
  CONTAINS
  
    SUBROUTINE CALCFORCE(NOPT,X,F,ETOT,ESCALE)
       USE MOD_BONDS, ONLY: ENERGY_BONDS
       USE MOD_ANGLES, ONLY: ENERGY_ANGLES
       USE MOD_DIHEDRALS, ONLY: ENERGY_DIHS
       USE MOD_DEBYEHUECKEL, ONLY: ENERGY_DH
       USE MOD_NONBONDED, ONLY: E_NONBONDED
       USE MOD_SAXS, ONLY: RNA_SAXS_FORCE
       USE MOD_RESTRAINTS, ONLY: E_DISTRESTR, E_POSRESTR, NRESTS, NPOSRES
       IMPLICIT NONE    
       INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
       REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
       REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
       REAL(KIND = REAL64), INTENT(OUT) :: ETOT
       REAL(KIND = REAL64), INTENT(IN), OPTIONAL :: ESCALE
       
       REAL(KIND = REAL64) :: THIS_E, THIS_F(NOPT), EHHB, ESTAK, EVDW, &
                              F_SAXS(NOPT), F_CONST(NOPT), ESAXS, ECONST
       
       ETOT = 0.0D0
       F(1:NOPT) = 0.0D0
       CALL RESET_POT_ENE(EVEC)
       !1. Bonded terms in potential
       CALL ENERGY_BONDS(NOPT, X, THIS_F, THIS_E)
       ETOT = ETOT + THIS_E
       F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)
       EVEC%EBOND = THIS_E
       !2. Bond angle terms
       CALL ENERGY_ANGLES(NOPT, X, THIS_F, THIS_E)
       ETOT = ETOT + THIS_E
       F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)
       EVEC%EANGLES = THIS_E 
       !3. Torsional energy
       CALL ENERGY_DIHS(NOPT, X, THIS_F, THIS_E)    
       ETOT = ETOT + THIS_E
       F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)
       EVEC%ETORS = THIS_E
       !4. Debye-Hueckel term
       CALL ENERGY_DH(NOPT, X, THIS_F, THIS_E)
       ETOT = ETOT + THIS_E
       F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)
       EVEC%EDH = THIS_E       
       !5.Non-bonded interactions
       CALL  E_NONBONDED(NOPT, X, THIS_F, EHHB, ESTAK, EVDW) 
       ETOT = ETOT + EHHB + ESTAK + EVDW
       F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)
       EVEC%EHBOND = EHHB
       EVEC%ESTAK = ESTAK  
       EVEC%EVDW = EVDW    
       !6. SAXS energy and force
       F_SAXS(1:NOPT) = 0.0D0
       CALL RNA_SAXS_FORCE(NOPT, X, ESAXS, F_SAXS)
       EVEC%ESAXS = ESAXS
       !7. Energy for any restraints
       ECONST = 0.0D0
       F_CONST(1:NOPT) = 0.0D0
       IF (NRESTS.GT.0) THEN
          CALL E_DISTRESTR(NOPT, X, THIS_F, THIS_E)
          ECONST = ECONST + THIS_E
          F_CONST(1:NOPT) = F_CONST(1:NOPT) + THIS_F(1:NOPT)
          EVEC%EDISTR = THIS_E
       ENDIF
       IF (NPOSRES.GT.0) THEN
          CALL E_POSRESTR(NOPT, X, THIS_F, THIS_E)
          ECONST = ECONST + THIS_E
          F_CONST(1:NOPT) = F_CONST(1:NOPT) + THIS_F(1:NOPT)
          EVEC%EPOSR = THIS_E
       ENDIF 
       ETOT = ETOT + ESAXS + ECONST
       F(1:NOPT) = F(1:NOPT) + F_SAXS(1:NOPT) + F_CONST(1:NOPT)
       EVEC%ETOT = ETOT                          
    END SUBROUTINE CALCFORCE 
    
    SUBROUTINE CALC_HBONDS(NOPT,X)
       USE HB_DEFS, ONLY: IHB, HBDAT, DO_HB, SAVE_HB
       USE MOD_NONBONDED, ONLY: E_NONBONDED
       IMPLICIT NONE    
       INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
       REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
       REAL(KIND = REAL64) :: EHB, ESTAK, EVDW, THIS_F(NOPT)
       
       IF (DO_HB) THEN
          IHB = IHB + 1
          SAVE_HB = .TRUE.
          WRITE(HBDAT, '(A,I6)') "#", IHB
          CALL E_NONBONDED(NOPT, X, THIS_F, EHB , ESTAK, EVDW)
          SAVE_HB = .FALSE.
       ELSE
         WRITE(*,*) "calc_forces> Error - Hydrogen-bond printing not initialised"
         STOP
       ENDIF
    END SUBROUTINE CALC_HBONDS

    !routine to reset energy vector
    SUBROUTINE RESET_POT_ENE(ENEPOT)
       TYPE(POT_ENE), INTENT(OUT) :: ENEPOT
       
       ENEPOT%EBOND = 0.0D0
       ENEPOT%EANGLES = 0.0D0 
       ENEPOT%EDH = 0.0D0       
       ENEPOT%ETORS = 0.0D0       
       ENEPOT%EHBOND = 0.0D0
       ENEPOT%ESTAK = 0.0D0  
       ENEPOT%EVDW = 0.0D0  
       ENEPOT%ESAXS = 0.0D0   
       ENEPOT%EDISTR = 0.0D0
       ENEPOT%EPOSR = 0.0D0
       ENEPOT%ETOT = 0.0D0 
    END SUBROUTINE RESET_POT_ENE

    !printing debug information
    SUBROUTINE PRINT_POT_ENE(EUNIT,ENEPOT)
       TYPE(POT_ENE), INTENT(IN) :: ENEPOT
       INTEGER, INTENT(IN) :: EUNIT
       
       WRITE(EUNIT, '(A,F15.5)') " Ebond:  ", ENEPOT%EBOND
       WRITE(EUNIT, '(A,F15.5)') " Eangle: ", ENEPOT%EANGLES
       WRITE(EUNIT, '(A,F15.5)') " Etors:  ", ENEPOT%ETORS
       WRITE(EUNIT, '(A,F15.5)') " Edh:    ", ENEPOT%EDH  
       WRITE(EUNIT, '(A,F15.5)') " Ehbond: ", ENEPOT%EHBOND
       WRITE(EUNIT, '(A,F15.5)') " Evdw:   ", ENEPOT%EVDW
       WRITE(EUNIT, '(A,F15.5)') " Estak:  ", ENEPOT%ESTAK
       WRITE(EUNIT, '(A,F15.5)') " Esaxs:  ", ENEPOT%ESAXS
       WRITE(EUNIT, '(A,F15.5)') " Edistr: ", ENEPOT%EDISTR
       WRITE(EUNIT, '(A,F15.5)') " Eposr:  ", ENEPOT%EPOSR
       WRITE(EUNIT, '(A,F15.5)') " Etot:   ", ENEPOT%ETOT         
    END SUBROUTINE PRINT_POT_ENE
    
END MODULE CALCFORCES

!TODO: need to add new routine/module to get Hbonds:

!         !lm759> compute and save Hbonds, only 
!         IF(do_hb.and.HBSAVET)THEN
!            save_hb=HBSAVET
!            ihb=ihb+1
!            write(hbdat,*) "#",ihb
!            write(7878,*) "#",ihb
!            call RNA_HYDROP(scale,X,F)
!            save_hb=.false.
!            return
!         ENDIF




