MODULE MOD_ANGLES
  USE PREC_HIRE
  INTEGER :: NANGLES  !Number of angles
  INTEGER :: NUMANG   !Number of unique angle types

  !Angles: E = TK * (T - TEQ)**2 with T the current angle
  REAL(KIND = REAL64), ALLOCATABLE :: TK(:)     !Angular force constants
  REAL(KIND = REAL64), ALLOCATABLE :: TEQ(:)    !Equilibrium angle  
  
  !book-keeping for angles
  INTEGER, ALLOCATABLE :: IT(:)                 !Angle atom1
  INTEGER, ALLOCATABLE :: JT(:)                 !Angle atom2
  INTEGER, ALLOCATABLE :: KT(:)                 !Angle atom3
  INTEGER, ALLOCATABLE :: ICT(:)                !Type of angle

  INTEGER, ALLOCATABLE :: NTHETATYPE(:)

  CONTAINS

   SUBROUTINE ALLOC_ANGLES()
      CALL DEALLOC_ANGLES()
      ALLOCATE(TK(NUMANG), TEQ(NUMANG), IT(NANGLES), JT(NANGLES), &
               KT(NANGLES), ICT(NANGLES), NTHETATYPE(NANGLES))
   END SUBROUTINE ALLOC_ANGLES  

   SUBROUTINE ASSIGN_THETATYPE()
      USE VAR_DEFS, ONLY: IAC
      IMPLICIT NONE
      INTEGER :: JN, I3, J3, K3T, L3T, IT0, IT1, IT2
      
      DO JN=1,NANGLES
         I3 = IT(JN)/3 + 1
         J3 = JT(JN)/3 + 1
         K3T = KT(JN)/3 + 1
         IT0 = IAC(I3)
         IT1 = IAC(J3)
         IT2 = IAC(K3T)

         IF (IT0.EQ.1) THEN
            !IT1=4 ; IT2=3
            IF (IT2.EQ.3) THEN
               NTHETATYPE(JN) = 4
            !IT1=4 ; IT2=5 
            ELSE IF (IT2.EQ.5) THEN
               NTHETATYPE(JN) = 6
            ENDIF
         ELSE IF (IT0.EQ.2) THEN
            !IT1=1 ; IT2=4
            IF (IT1.EQ.1) THEN
               NTHETATYPE(JN) = 3
            ENDIF        
         ELSE IF (IT0.EQ.3) THEN
            !IT1=2 ; IT2=1
            IF (IT1.EQ.2) THEN
               NTHETATYPE(JN) = 2
            ENDIF         
         ELSE IF (IT0.EQ.4) THEN
            !IT1=5 ; IT2=(6,8,10,11)
            IF (IT1.EQ.5) THEN
               NTHETATYPE(JN) = 0
            !IT1=3 ; IT2=2
            ELSE IF (IT1.EQ.3) THEN
               NTHETATYPE(JN) = 5              
            ENDIF         
         ELSE IF (IT0.EQ.5) THEN
            !IT1=6 ; IT2=7
            IF (IT1.EQ.6) THEN         
               NTHETATYPE(JN) = 1 
            !IT1=8 ; IT2=9 
            ELSE IF (IT1.EQ.8) THEN
               NTHETATYPE(JN) = 1 
            !IT1=4 ; IT2=3   
            ELSE IF (IT1.EQ.4) THEN
               NTHETATYPE(JN) = 7                    
            ENDIF
         ELSE
            WRITE(*,*) " assign_thetatype> ERROR: Unkown theta type"
            STOP
         ENDIF
      ENDDO
   END SUBROUTINE ASSIGN_THETATYPE


   !QUERY: Why is there an if TK(IC) .GE. 2.0D0?
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
      INTEGER :: P1, P2, P3

      !initialise force and energy
      F(1:NOPT) = 0.0d0
      EANGLE = 0.0d0

      DO JN = 1,NANGLES
         IC = ICT(JN)
         IF (TK(IC) .GE. 2.0d0) THEN
            I = IT(JN)
            J = JT(JN)
            K = KT(JN)
            RIJ = X(I+1:I+3)-X(J+1:J+3)
            RKJ = X(K+1:K+3)-X(J+1:J+3)

            RIJ0 = dot_product(RIJ, RIJ)
            RKJ0 = dot_product(RKJ, RKJ)
            RIK0 = dsqrt(RIJ0*RKJ0)
            CT0 = dot_product(RIJ, RKJ)/RIK0
            ! QUERY: This is clearly for numerical reasons, but it feels like 
            !        there ought to be a better solution for this using TANH?
            CT1 = MAX(-PT999,CT0)
            CT2 = MIN(PT999,CT1)
            ANT = DACOS(CT2)

            ! ENERGY
            DA = ANT-TEQ(IC)
            DF = TK(IC)*DA*SCORE_RNA(2+NTHETATYPE(JN))*SCORE_RNA(47)
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
         ENDIF
      END DO

   END SUBROUTINE ENERGY_ANGLES

   SUBROUTINE DEALLOC_ANGLES()
      IF (ALLOCATED(TK)) DEALLOCATE(TK)
      IF (ALLOCATED(TEQ)) DEALLOCATE(TEQ)
      IF (ALLOCATED(IT)) DEALLOCATE(IT)
      IF (ALLOCATED(JT)) DEALLOCATE(JT)
      IF (ALLOCATED(KT)) DEALLOCATE(KT)
      IF (ALLOCATED(ICT)) DEALLOCATE(ICT)
      IF (ALLOCATED(NTHETATYPE)) DEALLOCATE(NTHETATYPE)
   END SUBROUTINE DEALLOC_ANGLES

END MODULE MOD_ANGLES
