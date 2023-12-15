!> @file
!> Contains MOD_SAXS module handling the connection between the SAXS scoring routines and the energy calling interface

!> Module contains routine to obtain the SAXS force
MODULE MOD_SAXS
   USE PREC_HIRE
   USE SAXS_DEFS, ONLY: SAXSFORCET, SAXSSAVET, SAXSMODI, SAXSOFFI
   IMPLICIT NONE
 
   CONTAINS
      !> Subroutine to calculate SAXS force
      !> @brief
      !>
      !> Subroutine used to obtain SAXS data for a set of input coordinates.\n
      !> The routine calls the routines to generate and write the SAXS curves.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[out] ESAXS - SAXS energy
      !> @param[out] F - SAXS force
      !>
      !> @see GENERATE_SAXS_CURVE 
      !> @see WRITE_SAXS_CURVE_TO_UNIT
      SUBROUTINE RNA_SAXS_FORCE(NOPT, X, ESAXS, F, GRADT, HYDRATET)
         USE NAPARAMS, ONLY: SCORE_RNA
         USE NUM_DEFS, ONLY: PI
         USE SAXS_CALCS
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force
         REAL(KIND = REAL64), INTENT(OUT) :: ESAXS
         LOGICAL, INTENT(IN) :: GRADT !use SAXS force?
         LOGICAL, INTENT(IN) :: HYDRATET !use hydration routine

         REAL(KIND = REAL64) :: KSAXS, KDECRE, KMODUL
         REAL(KIND = REAL64), DIMENSION(NQPOINTS) :: LOGI

         ESAXS = 0.0D0
         F(1:NOPT) = 0.0D0
         CALC_SAXS_FORCE = COMPUTE_SAXS_SERIAL .AND. SAXSFORCET
         SAXS_SAVE = SAXS_PRINT .AND. SAXSSAVET
         
         IF (CALC_SAXS_FORCE) THEN
            CALC_SAXS_FORCE = .FALSE.
            KSAXS = SCORE_RNA(78)
            IF (MODULATE_SAXS_SERIAL) THEN
               !KMODUL = EXP(-(SIN(PI*SAXSOFFI)*SAXS_INVSIG)**2)
               KMODUL = 1
               !KDECRE = 0.5 * (1.0 + COS(2*PI*SAXSMODI))
               KDECRE = COS(PI*SAXSMODI)**4
            ELSE
               KMODUL = 1
               KDECRE = 1
            ENDIF
            ! WRITE(*,*) "KMODUL, KDECRE ", KMODUL, KDECRE
            !QUERY: Replace this magic number by a parameter or variable?
            IF (KMODUL * KDECRE .GE. 2.0D-2) THEN
               IF (HYDRATET) THEN
                  CALL GENERATE_SAXS_CURVE_HYDRATE(X, LOGI, ESAXS, F)
               ELSE
                  CALL GENERATE_SAXS_CURVE(X, LOGI, ESAXS, F, GRADT)
               END IF
               ESAXS = KSAXS*ESAXS*KMODUL*KDECRE
               F(1:NOPT) = KSAXS*KMODUL*KDECRE*F(1:NOPT)
            ELSE
               ESAXS = 0.0D0
            ENDIF
            IF (KMODUL * KDECRE .LE. 1.0D-10) THEN
               IF (HYDRATET) THEN
                  CALL GENERATE_SAXS_CURVE_HYDRATE(X, LOGI, ESAXS, F)
               ELSE
                  CALL GENERATE_SAXS_CURVE(X, LOGI, ESAXS, F, GRADT)
               END IF
               ESAXS = KSAXS*ESAXS
               CALL WRITE_SAXS_CURVE_TO_UNIT(LOGI,SAXSC)
               WRITE(*,*) 'SAXS score: ', ESAXS
               ESAXS = 0.0D0
            ENDIF
            IF (SAXS_PRINT .AND. KMODUL * KDECRE .GE. 2.0D-2) THEN
               CALL WRITE_SAXS_CURVE_TO_UNIT(LOGI,SAXSC)
               !WRITE(SAXSS, *) KMODUL*KDECRE, ESAXS
               WRITE(*,*) 'SAXS score: ', ESAXS, KMODUL*KDECRE
            END IF
            !WRITE(*,*) 'SAXS score: ', ESAXS
            !WRITE(*,*) 'SAXS force: ', F
         END IF 
      END SUBROUTINE RNA_SAXS_FORCE

END MODULE MOD_SAXS
