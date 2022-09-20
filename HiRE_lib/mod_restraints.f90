!> @file
!> Contains MOD_RESTRAINTS to handle distance and position restraints

!> Module containing the required variables and functions to deal with restraints
MODULE MOD_RESTRAINTS
   USE PREC_HIRE
   ! Variables for restraining RNA
   !> Number of distance restraints
   INTEGER :: NRESTS                         
   !> Number of position restraints
   INTEGER :: NPOSRES                        
   ! Distance restraints
   !> Restraint atom1 (distance restraint)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: RESTI
   !> Restraint atom2 (distance restraint)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: RESTJ 
   !> Restraint type (distance restraint)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: TREST 

   !> Restraint force constant (distance restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: RESTK 
   !> Length (distance restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: RESTL   
   !> Restraint step to update restraints (distance restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: DRESTL 

   ! Position restraints
   !> Constraint ID (position restraint)
   INTEGER, ALLOCATABLE, DIMENSION(:) :: PRI 
   !> Constraint force (position restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: PRK
   !> Restraint coordinates (position restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:,:) :: PRX  
   !> Values to update restraint (position restraint)
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:,:) :: PRDX 
   
   CONTAINS
      !> Allocate distance restraints
      SUBROUTINE ALLOC_DISTRESTR()
         CALL DEALLOC_DISTRESTR()
         ALLOCATE(RESTI(NRESTS), RESTJ(NRESTS), TREST(NRESTS), &
                  RESTK(NRESTS), RESTL(NRESTS), DRESTL(NRESTS))
      END SUBROUTINE ALLOC_DISTRESTR

      !> Allocate position restraint
      SUBROUTINE ALLOC_POSRESTR()
         CALL DEALLOC_POSRESTR()
         ALLOCATE(PRI(NPOSRES), PRK(NPOSRES), PRX(3,NPOSRES), PRDX(3,NPOSRES))
      END SUBROUTINE ALLOC_POSRESTR

      !> Calculate energy for distance restraints
      !> @brief
      !>
      !> Routine to calculate the distance restraints' contributions to the energy and force.\n
      !> The restraints can be modulated through steptime and updated through DRESTL.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[out] F - forces from restraint
      !> @param[out] EREST - restraint energy
      !> @param[in] STEPTIME - modulation of restraints (COS(STEPTIME/PI/5)**2)
      SUBROUTINE E_DISTRESTR(NOPT, X, F, EREST, STEPTIME)
         USE NUM_DEFS, ONLY: PI
         USE VEC_UTILS, ONLY: EUC_NORM
         IMPLICIT NONE
      
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EREST
         REAL(KIND = REAL64), INTENT(IN), OPTIONAL :: STEPTIME

         REAL(KIND = REAL64) :: ECON, ESCALE, DIFF(3), DLEN, V, K, DF(3)
         INTEGER :: IDX, I, J
   
         IF (PRESENT(STEPTIME)) THEN
            ESCALE = COS(STEPTIME/PI/5)**2
         ELSE
            ESCALE = 1.0D0
         ENDIF
         
         EREST = 0.0D0
      
         DO IDX = 1, NRESTS
            I = RESTI(IDX)
            J = RESTJ(IDX)
         
            DIFF(1:3) = X((3*I-2):(3*I)) - X((3*J-2):(3*J))
            DLEN = EUC_NORM(DIFF)
         
            ! Linear for d>2.0
            V = DLEN - RESTL(IDX)
            K = RESTK(IDX)
            IF ((TREST(IDX).GE.2).AND.(ABS(V).GT.2.0D0)) THEN
               ECON = 4.0D0*K*(ABS(V) - 1.0D0)
               DF(1:3) = 4.0D0*K*V/(ABS(V)*DLEN)*DIFF(1:3)
            ELSE 
               ECON = K*V**2
               DF(1:3) = 2*K*V/DLEN*DIFF(1:3)
            ENDIF
            IF (TREST(IDX).EQ.2) THEN
               ECON = ECON*ESCALE
               DF(1:3) = ESCALE*DIFF(1:3)
            ENDIF
            !add contribution to E and F
            EREST = EREST + ECON
            F((3*I-2):(3*I)) = F((3*I-2):(3*I)) - DF(1:3)
            F((3*J-2):(3*J)) = F((3*J-2):(3*J)) + DF(1:3)
            !update restraint
            RESTL(IDX) = RESTL(IDX) + DRESTL(IDX)
         ENDDO
      END SUBROUTINE E_DISTRESTR
   
      !> Calculate energy for position restraints
      !> @brief
      !>
      !> Routine to calculate the position restraints' contributions to the energy and force.\n
      !> The restraints can be updated through PRDX.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[out] F - forces from restraint
      !> @param[out] EREST - restraint energy
      SUBROUTINE E_POSRESTR(NOPT, X, F, EREST)
         IMPLICIT NONE     
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EREST
   
         REAL(KIND = REAL64) :: DIFF(3), V, K, ECON, DF(3)
         INTEGER :: IDX, I
         
         EREST = 0.0D0
         
         DO IDX = 1,NPOSRES
            I = PRI(IDX)
            DIFF(1:3) = X((3*I-2):(3*I)) - PRX(1:3,IDX)
            V = DOT_PRODUCT(DIFF,DIFF)
            K = PRK(IDX)
            EREST = EREST + K*V
            DF(1:3) = 2.0D0*K*DIFF(1:3)
            F((3*I-2):(3*I)) = F((3*I-2):(3*I)) - DF(1:3)
         ENDDO
         IF (NPOSRES.GT.0) THEN
            PRX = PRX + PRDX
         ENDIF  
      END SUBROUTINE E_POSRESTR

      !> Deallocate distance restraint arrays
      SUBROUTINE DEALLOC_DISTRESTR()
         IF (ALLOCATED(RESTI)) DEALLOCATE(RESTI)
         IF (ALLOCATED(RESTJ)) DEALLOCATE(RESTJ)
         IF (ALLOCATED(TREST)) DEALLOCATE(TREST)
         IF (ALLOCATED(RESTK)) DEALLOCATE(RESTK)
         IF (ALLOCATED(RESTL)) DEALLOCATE(RESTL)
         IF (ALLOCATED(DRESTL)) DEALLOCATE(DRESTL)
      END SUBROUTINE DEALLOC_DISTRESTR

      !> Deallocate position restraint arrays
      SUBROUTINE DEALLOC_POSRESTR()
         IF (ALLOCATED(PRI)) DEALLOCATE(PRI)
         IF (ALLOCATED(PRK)) DEALLOCATE(PRK)
         IF (ALLOCATED(PRX)) DEALLOCATE(PRX)
         IF (ALLOCATED(PRDX)) DEALLOCATE(PRDX)
      END SUBROUTINE DEALLOC_POSRESTR
   
END MODULE MOD_RESTRAINTS

