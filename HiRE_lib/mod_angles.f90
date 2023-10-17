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

   !Angles: E = TK * (T - TEQ)**2 with T the current angle
   !> Angular force constants
   REAL(KIND = REAL64), ALLOCATABLE :: TK(:)     
   !> Equilibrium angle  
   REAL(KIND = REAL64), ALLOCATABLE :: TEQ(:)  
  
   !book-keeping for angles
   !> Angle atom1
   INTEGER, ALLOCATABLE :: IT(:)                 
   !> Angle atom2
   INTEGER, ALLOCATABLE :: JT(:)                 
   !> Angle atom3
   INTEGER, ALLOCATABLE :: KT(:)                 
   !> Type of angle
   INTEGER, ALLOCATABLE :: ICT(:)                

   !> Angle type information needed for CG scaling
   INTEGER, ALLOCATABLE :: NTHETATYPE(:)

   CONTAINS
      !> Routine to allocate all required arrays
      SUBROUTINE ALLOC_ANGLES()
         CALL DEALLOC_ANGLES()
         ALLOCATE(TK(NUMANG), TEQ(NUMANG), IT(NANGLES), JT(NANGLES), &
                  KT(NANGLES), ICT(NANGLES), NTHETATYPE(NANGLES))
      END SUBROUTINE ALLOC_ANGLES  

      !> Assign the angle type to obtain the correct CG scaling information
      SUBROUTINE ASSIGN_THETATYPE()
         USE VAR_DEFS, ONLY: IAC, RESTYPES, ATOMINRES_LOOKUP
         IMPLICIT NONE
         INTEGER :: JN, I3, J3, K3T, IT0, IT1, IT2, MOLTYPE

         ! Atom types (IAC):
         ! 1: C  2: O  3: P  4: R4  5: R1  6: S1
         ! 7,8: G1,2  9,10: A1,2  11: U1  12: C1
         ! 13: D 14: MG 15: NA  16:CL
         ! RNA angle types (Theta type, angle, IACs)
         ! 0 - R4-R1-X1         4-5-(7,9,11,12)
         ! 1 - R1-A1/G1-A2/G2   5-(7/9)-(8/10)
         ! 2 - P-O-C            3-2-1
         ! 3 - O-C-R4           2-1-4
         ! 4 - C-R4-P           1-4-3
         ! 5 - R4-P-O           4-3-2
         ! 6 - C-R4-R1          1-4-5
         ! 7 - R1-R4-P          5-4-3
         ! DNA angle types (Theta type, angle, IACs)
         ! 8 - R4-S1-X1         4-6-(7,9,11,12)
         ! 9 - S1-A1/G1-A2/G2   6-(7/9)-(8/10)
         !10 - P-O-C            3-2-1
         !11 - O-C-R4           2-1-4
         !12 - C-R4-P           1-4-3
         !13 - R4-P-O           4-3-2
         !14 - C-R4-S1          1-4-6
         !15 - S1-R4-P          6-4-3    
         !
         !QUERY: Potentially we could merge 2 and 10, 3 and 11, 4 and 12 and 5 and 13
         !       Will they actually be different?     
         DO JN=1,NANGLES
            I3 = IT(JN)/3 + 1
            J3 = JT(JN)/3 + 1
            K3T = KT(JN)/3 + 1
            IT0 = IAC(I3)
            IT1 = IAC(J3)
            IT2 = IAC(K3T)
            ! MOLTYPE is coding for RNA (0) or DNA (1)
            MOLTYPE = RESTYPES(ATOMINRES_LOOKUP(K3T))

            IF (IT0.EQ.1) THEN
               !IT1=4 ; IT2=3
               IF (IT2.EQ.3) THEN
                  IF (MOLTYPE.EQ.0) THEN
                     NTHETATYPE(JN) = 4
                  ELSE IF (MOLTYPE.EQ.1) THEN
                     NTHETATYPE(JN) = 12
                  END IF
               !IT1=4 ; IT2=5 
               ELSE IF (IT2.EQ.5) THEN                  
                  NTHETATYPE(JN) = 6
               !IT1=4 ; IT2=6  
               ELSE IF (IT2.EQ.6) THEN                  
                  NTHETATYPE(JN) = 14
               ENDIF
            ELSE IF (IT0.EQ.2) THEN
               !IT1=1 ; IT2=4
               IF (IT1.EQ.1) THEN
                  IF (MOLTYPE.EQ.0) THEN
                     NTHETATYPE(JN) = 3
                  ELSE IF (MOLTYPE.EQ.1) THEN
                     NTHETATYPE(JN) = 11
                  END IF
               ENDIF        
            ELSE IF (IT0.EQ.3) THEN
               !IT1=2 ; IT2=1
               IF (IT1.EQ.2) THEN
                  IF (MOLTYPE.EQ.0) THEN
                     NTHETATYPE(JN) = 2
                  ELSE IF (MOLTYPE.EQ.1) THEN
                     NTHETATYPE(JN) = 10
                  END IF                  
               ENDIF         
            ELSE IF (IT0.EQ.4) THEN
               !IT1=5 ; IT2=(6,8,10,11)
               IF (IT1.EQ.5) THEN
                  NTHETATYPE(JN) = 0
               ELSE IF (IT1.EQ.6) THEN
                  NTHETATYPE(JN) = 8
               !IT1=3 ; IT2=2
               ELSE IF (IT1.EQ.3) THEN
                  IF (MOLTYPE.EQ.0) THEN
                     NTHETATYPE(JN) = 5
                  ELSE IF (MOLTYPE.EQ.1) THEN
                     NTHETATYPE(JN) = 13
                  END IF           
               ENDIF         
            ELSE IF (IT0.EQ.5) THEN
               !IT1=7 ; IT2=8
               IF (IT1.EQ.7) THEN         
                  NTHETATYPE(JN) = 1 
               !IT1=9 ; IT2=10 
               ELSE IF (IT1.EQ.9) THEN
                  NTHETATYPE(JN) = 1 
               !IT1=4 ; IT2=3   
               ELSE IF (IT1.EQ.4) THEN
                  NTHETATYPE(JN) = 7                    
               ENDIF
            ELSE IF (IT0.EQ.6) THEN
               !IT1=7 ; IT2=8
               IF (IT1.EQ.7) THEN         
                  NTHETATYPE(JN) = 9 
               !IT1=9 ; IT2=10 
               ELSE IF (IT1.EQ.9) THEN
                  NTHETATYPE(JN) = 9 
               !IT1=4 ; IT2=3   
               ELSE IF (IT1.EQ.4) THEN
                  NTHETATYPE(JN) = 15                    
               ENDIF               
            ELSE
               WRITE(*,*) " assign_thetatype> ERROR: Unkown theta type"
               STOP
            ENDIF
         ENDDO
      END SUBROUTINE ASSIGN_THETATYPE

      !> Calculate the energy and force contribution from the angular terms
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
               DF = TK(IC)*DA*SCORE_RNA(3+NTHETATYPE(JN))*SCORE_RNA(2)
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

      !> Deallocate all arrays in this module
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
