!> @file
!> Contains MOD_ANGLES to handle bonding angles

!> Module containing all routines and variables to calculate the angular contribution to E and F\n
!> Angles: E = TK * (T - TEQ)**2 with T the current angle
MODULE MOD_ANGLES
   USE PREC_HIRE
   !> Number of angles
   INTEGER :: NANGLES  
   !> Number of unique angle types
   INTEGER :: NUMANG   

   !> Number of angles
   INTEGER :: NQANGLES  
   !> Number of unique angle types
   INTEGER :: NUMQANG   

   !Angles: E = TK * (T - TEQ)**2 with T the current angle
   !> Angular force constants
   REAL(KIND = REAL64), ALLOCATABLE :: TK(:)     
   !> Equilibrium angle  
   REAL(KIND = REAL64), ALLOCATABLE :: TEQ(:)  
  
   !book-keeping for angles - square potential (one equilibrium angle)
   !> Angle atom1
   INTEGER, ALLOCATABLE :: IT(:)                 
   !> Angle atom2
   INTEGER, ALLOCATABLE :: JT(:)                 
   !> Angle atom3
   INTEGER, ALLOCATABLE :: KT(:)                 
   !> Type of angle
   INTEGER, ALLOCATABLE :: ICT(:)  
   
   !Quartic angles
   ! "force" constant of overall energy (effectively scaling)
   REAL(KIND = REAL64), ALLOCATABLE :: TKQ(:)
   ! Equilibrium angles (correpsond to the maximum in the quartic potential)
   REAL(KIND = REAL64), ALLOCATABLE :: REFQ(:)

   !book keeping
   !> Angle atom 1
   INTEGER, ALLOCATABLE :: ITQ(:)                 
   !> Angle atom2
   INTEGER, ALLOCATABLE :: JTQ(:)                 
   !> Angle atom3
   INTEGER, ALLOCATABLE :: KTQ(:)  
   !> Type of qangle
   INTEGER, ALLOCATABLE :: ICTQ(:)  

   !potential parameters
   REAL(KIND = REAL64), ALLOCATABLE :: AQ(:), BQ(:), CQ(:), DQ(:), EQ(:)

   CONTAINS
      !> Routine to allocate all required arrays
      SUBROUTINE ALLOC_ANGLES()
         CALL DEALLOC_ANGLES()
         ALLOCATE(TK(NUMANG), TEQ(NUMANG), IT(NANGLES), JT(NANGLES), &
                  KT(NANGLES), ICT(NANGLES))
         ALLOCATE(TKQ(NUMQANG), REFQ(NUMQANG), AQ(NUMQANG), BQ(NUMQANG), &
                  CQ(NUMQANG), DQ(NUMQANG), EQ(NUMQANG), ITQ(NQANGLES), &
                  JTQ(NQANGLES), KTQ(NQANGLES), ICTQ(NQANGLES))
      END SUBROUTINE ALLOC_ANGLES  

      !> Routine to get all anglar contributions
      SUBROUTINE ENERGY_ALL_ANGLES(NOPT, X, F, EANGLE)
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EANGLE  

         REAL(KIND = REAL64) :: F1(NOPT), F2(NOPT), E1, E2

         !initialise force and energy
         F(1:NOPT) = 0.0d0
         EANGLE = 0.0d0

         CALL ENERGY_ANGLES(NOPT, X, F1, E1)
         CALL ENERGY_QANGLES(NOPT, X, F2, E1)

         EANGLE = E1 + E2
         F(1:NOPT) = F1(1:NOPT) + F2(1:NOPT)
      END SUBROUTINE ENERGY_ALL_ANGLES

      !> Calculate the energy and force contribution from the angular terms
      SUBROUTINE ENERGY_ANGLES(NOPT, X, F, EANGLE)
         USE NAPARAMS, ONLY: SCORE_RNA
         USE PREC_HIRE
         IMPLICIT NONE  
      
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EANGLE  

         REAL(KIND = REAL64) :: RIJ(3), RKJ(3), RIJ0, RKJ0, RIK0
         REAL(KIND = REAL64) :: RDI(3), RDJ(3), RDK(3)
         REAL(KIND = REAL64) :: CT0, CT1, CT2, ANT, DA, DF, EAW, DFW
         REAL(KIND = REAL64), PARAMETER :: PT999 = 0.999d0
         INTEGER :: JN, I, J, K, IC

         !initialise force and energy
         F(1:NOPT) = 0.0d0
         EANGLE = 0.0d0

         DO JN = 1,NANGLES
            IC = ICT(JN)
            I = IT(JN)
            J = JT(JN)
            K = KT(JN)
            RIJ = X(I+1:I+3)-X(J+1:J+3)
            RKJ = X(K+1:K+3)-X(J+1:J+3)

            RIJ0 = dot_product(RIJ, RIJ)
            RKJ0 = dot_product(RKJ, RKJ)
            RIK0 = dsqrt(RIJ0*RKJ0)
            CT0 = dot_product(RIJ, RKJ)/RIK0
            CT1 = MAX(-PT999,CT0)
            CT2 = MIN(PT999,CT1)
            ANT = DACOS(CT2)

            ! ENERGY
            DA = ANT-TEQ(IC)
            DF = TK(IC)*DA*SCORE_RNA(2)
            EAW = DF*DA
            DFW = -(2*DF)/DSIN(ANT)

            EANGLE = EANGLE + EAW
            ! FORCE
            rDI = DFW*(rKJ/RIK0-CT2*rIJ/RIJ0)
            rDK = DFW*(rIJ/RIK0-CT2*rKJ/RKJ0)
            rDJ = -rDI-rDK
            F(I+1:I+3) = F(I+1:I+3) - rDI
            F(J+1:J+3) = F(J+1:J+3) - rDJ
            F(K+1:K+3) = F(K+1:K+3) - rDK                 
         END DO

      END SUBROUTINE ENERGY_ANGLES

      SUBROUTINE ENERGY_QANGLES(NOPT, X, F, EQANGLE)
         USE NAPARAMS, ONLY: SCORE_RNA
         USE PREC_HIRE
         IMPLICIT NONE  
      
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EQANGLE  

         REAL(KIND = REAL64) :: RIJ(3), RKJ(3), RIJ0, RKJ0, RIK0
         REAL(KIND = REAL64) :: RDI(3), RDJ(3), RDK(3)
         REAL(KIND = REAL64) :: CT0, CT1, CT2, ANT
         REAL(KIND = REAL64) :: A, B, C, D, E, REF, KANG, DIFF, DIFF2, DIFF3, DIFF4, EQDUMMY, DF
         REAL(KIND = REAL64), PARAMETER :: PT999 = 0.999d0
         INTEGER :: JN, I, J, K, IC

         !initialise force and energy
         F(1:NOPT) = 0.0d0
         EQANGLE = 0.0d0

         DO JN = 1,NQANGLES
            I = ITQ(JN)
            J = JTQ(JN)
            K = KTQ(JN)
            IC = ICTQ(JN)
            RIJ = X(I+1:I+3)-X(J+1:J+3)
            RKJ = X(K+1:K+3)-X(J+1:J+3)

            RIJ0 = dot_product(RIJ, RIJ)
            RKJ0 = dot_product(RKJ, RKJ)
            RIK0 = dsqrt(RIJ0*RKJ0)
            CT0 = dot_product(RIJ, RKJ)/RIK0
            CT1 = MAX(-PT999,CT0)
            CT2 = MIN(PT999,CT1)
            ANT = DACOS(CT2)

            !Energy of angle term
            A = AQ(IC)
            B = BQ(IC)
            C = CQ(IC)
            D = DQ(IC)
            E = EQ(IC)
            REF = REFQ(IC)
            KANG = TKQ(IC)

            DIFF = ANT - REF
            DIFF2 = DIFF**2
            DIFF3 = DIFF**3
            DIFF4 = DIFF**4

            EQDUMMY = KANG*(A*DIFF4 + B*DIFF3 + C*DIFF2 + D*DIFF + E)*SCORE_RNA(2)
            EQANGLE = EQANGLE + EQDUMMY

            !Force
            DF = KANG*(4*A*DIFF3+3*B*DIFF2 + 2*C*DIFF + D)*SCORE_RNA(2) ! derivative of V_qangle with respect to theta
            ! now need derivative of theta with repsect to x,y,z
            rDI = DF*(rKJ/RIK0-CT2*rIJ/RIJ0)
            rDK = DF*(rIJ/RIK0-CT2*rKJ/RKJ0)
            rDJ = -rDI-rDK
            F(I+1:I+3) = F(I+1:I+3) - rDI
            F(J+1:J+3) = F(J+1:J+3) - rDJ
            F(K+1:K+3) = F(K+1:K+3) - rDK  
         END DO
      END SUBROUTINE ENERGY_QANGLES

      !> Deallocate all arrays in this module
      SUBROUTINE DEALLOC_ANGLES()
         IF (ALLOCATED(TK)) DEALLOCATE(TK)
         IF (ALLOCATED(TEQ)) DEALLOCATE(TEQ)
         IF (ALLOCATED(IT)) DEALLOCATE(IT)
         IF (ALLOCATED(JT)) DEALLOCATE(JT)
         IF (ALLOCATED(KT)) DEALLOCATE(KT)
         IF (ALLOCATED(ICT)) DEALLOCATE(ICT)
         IF (ALLOCATED(TKQ)) DEALLOCATE(TKQ)
         IF (ALLOCATED(REFQ)) DEALLOCATE(REFQ)
         IF (ALLOCATED(AQ)) DEALLOCATE(AQ)
         IF (ALLOCATED(BQ)) DEALLOCATE(BQ)
         IF (ALLOCATED(CQ)) DEALLOCATE(CQ)
         IF (ALLOCATED(DQ)) DEALLOCATE(DQ)
         IF (ALLOCATED(EQ)) DEALLOCATE(EQ)
         IF (ALLOCATED(ITQ)) DEALLOCATE(ITQ)
         IF (ALLOCATED(JTQ)) DEALLOCATE(JTQ)
         IF (ALLOCATED(KTQ)) DEALLOCATE(KTQ)
         IF (ALLOCATED(ICTQ)) DEALLOCATE(ICTQ)
      END SUBROUTINE DEALLOC_ANGLES

END MODULE MOD_ANGLES
