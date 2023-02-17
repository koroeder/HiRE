!> @file 
!> Contains CALCFORCES module

!> Module to calculate forces and energies
!> @brief
!>
!> Module containsing the subroutines that interface the individual energy components.
MODULE CALCFORCES
   USE PREC_HIRE
   IMPLICIT NONE

   !> Potential energy components
   TYPE POT_ENE
      !> bonding energy
      REAL(KIND = REAL64) :: EBOND     
      !> bond angles
      REAL(KIND = REAL64) :: EANGLES
      !> torsional energy
      REAL(KIND = REAL64) :: ETORS         
      !> Debye-Hueckel energy 
      REAL(KIND = REAL64) :: EDH           
      !> hydrogen bonding
      REAL(KIND = REAL64) :: EHBOND        
      !> excluded volume
      REAL(KIND = REAL64) :: EVDW           
      !> stacking interactions
      REAL(KIND = REAL64) :: ESTAK         
      !> SAXS energy
      REAL(KIND = REAL64) :: ESAXS          
      !> Distance contraints
      REAL(KIND = REAL64) :: EDISTR 
      !> Positional contraints
      REAL(KIND = REAL64) :: EPOSR          
      !> Total energy 
      REAL(KIND = REAL64) :: ETOT          
   END TYPE POT_ENE

   !> vector containing the energy contributions
   TYPE(POT_ENE) :: EVEC

   !> Individual energy scaling - can be used for H-REX simulation focusing for example on non-bonding terms
   REAL(KIND=REAL64) :: SCALING(7) = 1.0D0

   CONTAINS
  
      !> Subroutine calculating the total energy and gradient
      !> @brief
      !>
      !> The subroutine first clears the energy values saved by resetting EVEC.\n
      !> Then, in turn, all the energy components are called.\n
      !> After each call the new energy component is saved, and energy and gradient added to the total values.
      !>
      !> @param[in] NOPT - number of variables (three times the number of particles)
      !> @param[in] X - coordinates
      !> @param[out] ETOT - total energy
      !> @param[out] F - gradient
      !> @param[in] ESCALE - optional scaling parameter, currently not used
      SUBROUTINE CALCFORCE(NOPT,X,F,ETOT,ESCALE)
         USE MOD_BONDS, ONLY: ENERGY_BONDS
         USE MOD_ANGLES, ONLY: ENERGY_ANGLES
         USE MOD_DIHEDRALS, ONLY: ENERGY_DIHS
         USE MOD_DEBYEHUECKEL, ONLY: ENERGY_DH, DH_ENERGY
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
         ETOT = ETOT + THIS_E*SCALING(1)
         F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)*SCALING(1)
         EVEC%EBOND = THIS_E*SCALING(1)
         !2. Bond angle terms
         CALL ENERGY_ANGLES(NOPT, X, THIS_F, THIS_E)
         ETOT = ETOT + THIS_E*SCALING(2)
         F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)*SCALING(2)
         EVEC%EANGLES = THIS_E*SCALING(2) 
         !3. Torsional energy
         CALL ENERGY_DIHS(NOPT, X, THIS_F, THIS_E)    
         ETOT = ETOT + THIS_E*SCALING(3)
         F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)*SCALING(3)
         EVEC%ETORS = THIS_E*SCALING(3)
         !4. Debye-Hueckel term
         !CALL ENERGY_DH(NOPT, X, THIS_F, THIS_E)
         CALL DH_ENERGY(NOPT, X, THIS_F, THIS_E)          
         ETOT = ETOT + THIS_E*SCALING(4)
         F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)*SCALING(4)
         EVEC%EDH = THIS_E*SCALING(4)  
         WRITE(*,'(A,F20.10,A,F20.10)') "DH E: ", THIS_E, "  Force: ", MAX(SQRT(SUM(THIS_F(1:NOPT)**2)/NOPT), 1.0D-100 )
         !5.Non-bonded interactions
         CALL  E_NONBONDED(NOPT, X, THIS_F, EHHB, ESTAK, EVDW) 
         ETOT = ETOT + EHHB*SCALING(5) + ESTAK*SCALING(6) + EVDW*SCALING(7)
         ! not ideal but I guess the average is the best with the current setup, probably should make the scaling identical for all three
         F(1:NOPT) = F(1:NOPT) + THIS_F(1:NOPT)*(SCALING(5)+SCALING(6)+SCALING(7))/3.0
         EVEC%EHBOND = EHHB*SCALING(5)
         EVEC%ESTAK = ESTAK*SCALING(6) 
         EVEC%EVDW = EVDW*SCALING(7)   
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
      
      !> Get hydrogen-bonding information
      !> @brief
      !>
      !> If hydrogen bonding has been initialised, this calls the hydrogen-bonding energies 
      !> through the nonbonded energy contribution.
      !> A file with the hydrogen bonding information is returned.
      !>
      !> @param[in] NOPT - number of variables
      !> @param[in] X - coordinates
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

      !> Routine to reset energy vector
      !> @brief
      !>
      !> Resets all energy components of the provided energy vector.
      !>
      !> @param[in] ENEPOT - energy vector to be reset
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

      !> Printing debug information
      !> @brief
      !>
      !> Subroutines prints all energy components to an output file.
      !>
      !> @param[in] EUNIT - output file unit
      !> @param[in] ENEPOT - energy vector to be printed out
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




