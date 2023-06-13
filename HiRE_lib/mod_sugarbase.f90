!> @file
!> Contains MOD_SUGARBASE handling the sugar base interaction term

!> Module containing all routines and variables to calculate the sugar-base contributions
MODULE MOD_SUGARBASE
   USE PREC_HIRE

   !> Number of SB pairs
   INTEGER :: NSBPAIRS = 0
   !> List of pairs considered 
   INTEGER, ALLOCATABLE :: SB_PAIRS(:,:)
   !> scaling for sugar base interactions
   REAL(KIND = REAL64) :: SB_SCALE
   !> potential parameters for LJ potential
   REAL(KIND = REAL64) :: ALJ12
   REAL(KIND = REAL64) :: BLJ6

   CONTAINS

      SUBROUTINE INIT_SUGARBASE()
         USE VAR_DEFS, ONLY: NRES, RESTYPES, RESSTART, RESFINAL, IGRAPH
         USE NAPARAMS, ONLY: SCORE_RNA
         IMPLICIT NONE
         INTEGER :: DUMMY_PAIRS(NRES,2)
         INTEGER :: I, J, FIRST, LAST, IDBASE, IDSUGAR
         CHARACTER(4) :: ATNAME

         SB_SCALE = SCORE_RNA(55)
         ALJ12 = SCORE_RNA(56)
         BLJ6 = SCORE_RNA(57)

         CALL DEALLOC_SB()

         DO I=1,NRES
            !only for RNA
            IF (RESTYPES(I).EQ.0) THEN
               FIRST = RESSTART(I)
               LAST = RESFINAL(I)
               DO J=FIRST,LAST
                  ATNAME = IGRAPH(J)
                  IF (ATNAME.EQ."R4") IDSUGAR = J
                  IF ((ATNAME.EQ."A1").OR.(ATNAME.EQ."C1").OR.(ATNAME.EQ."G1").OR. &
                      (ATNAME.EQ."U1")) THEN
                     IDBASE = J
                  END IF
               END DO
               NSBPAIRS = NSBPAIRS + 1
               DUMMY_PAIRS(NSBPAIRS,1) = IDSUGAR
               DUMMY_PAIRS(NSBPAIRS,2) = IDBASE
            END IF
         END DO
         ALLOCATE(SB_PAIRS(NSBPAIRS,2))
         SB_PAIRS(1:NSBPAIRS,1:2) = DUMMY_PAIRS(1:NSBPAIRS,1:2)
      END SUBROUTINE INIT_SUGARBASE

      SUBROUTINE E_SUGARBASE(NOPT, X, F, ESB)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: ESB
         INTEGER :: I, J, AT1, AT2
         REAL(KIND = REAL64) :: DIST2, DIST6, DIST12, DIST8, DIST14, FA, FSB(3)

         ESB = 0.0
         F(1:NOPT) = 0.0D0
         DO I=1,NSBPAIRS
            AT1 = SB_PAIRS(I,1)
            AT2 = SB_PAIRS(I,2)
            !energy
            DIST2 = (X(3*AT1-2)-X(3*AT2-2))**2 + (X(3*AT1-1)-X(3*AT2-1))**2 + (X(3*AT1)-X(3*AT2))**2
            DIST6 = DIST2**3
            DIST12 = DIST6**2
            ESB = ESB + ALJ12/DIST12 - BLJ6/DIST6
            !forces
            DIST8 = DIST2*DIST6
            DIST14 = DIST8*DIST6
            FA = 12.0*ALJ12/DIST14 - 6.0*BLJ6/DIST8
            DO J=1,3
               FSB(J) = FA * (X(3*(AT1-1)+J)-X(3*(AT2-1)+J)) 
            END DO
            F(3*AT1-2:3*AT1) = F(3*AT1-2:3*AT1) + FSB(1:3)
            F(3*AT2-2:3*AT2) = F(3*AT2-2:3*AT2) - FSB(1:3)
         END DO
         ESB = ESB * SB_SCALE
         F(1:NOPT) = F(1:NOPT) * SB_SCALE
      END SUBROUTINE E_SUGARBASE

      SUBROUTINE DEALLOC_SB()
         IF (ALLOCATED(SB_PAIRS)) DEALLOCATE(SB_PAIRS)
      END SUBROUTINE DEALLOC_SB
END MODULE MOD_SUGARBASE