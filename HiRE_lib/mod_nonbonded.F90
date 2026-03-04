!> @file
!> Contains MOD_NONBONDED to calculate non-bonded contributions

!> Module to calculate all non-bonded contributions including the hydrogen-bonding, excluded volume and stacking energy
MODULE MOD_NONBONDED
   USE PREC_HIRE
   USE NBDEFS
   
   CONTAINS
      !> Routine to calaculate non-bonded energy
      !> @brief
      !>
      !> This subroutine calculates the non-bonded terms for all residues.\n
      !> An iteration of all base pairs takes place. Firstly, if the bases are far apart, no further calculations are conducted.\n
      !> Currently, the cutoff value for this is 20 Angstrom and hard-coded.\n
      !> If the bases are close enough, the hydrogen bonding and stacking energies are computed first.\n
      !> Then the excluded volume between particles is calaculated.
      !>
      !> @warning The cutoff for calculations is hard coded here to the distance squared being less than 400.
      !>
      !> @param[in] NOPT - number of degrees of freedom
      !> @param[in] X - input coordinates
      !> @param[out] F - gradient array for the non-bonded contributions
      !> @param[out] EHHB - hydrogen bonding energy
      !> @param[out] ESTAK - base stacking energy
      !> @param[out] EVDW - excluded volume interactions 
      !>
      !> @see RNA_BB
      !> @see RNA_STACKV2
      !> @see ENERGY_EXCLV
      SUBROUTINE E_NONBONDED(NOPT, X, F, EHHB, ESTAK, EVDW)
         USE UTILS_IO, ONLY: GETUNIT
         USE MOD_HBONDS, ONLY: ENERGY_HB
         USE MOD_EXCLV, ONLY: ENERGY_EXV
         USE MOD_BASESTACKING, ONLY: NA_STACKV2
         USE NAPARAMS, ONLY: BOCC, BTYPE, BLIST, RCUT2_EXCLV, NBCUT
         USE VAR_DEFS, ONLY: NRES, RESSTART, RESFINAL, RESTYPES, IAC
         IMPLICIT NONE
      
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force 
         REAL(KIND = REAL64), INTENT(OUT) :: EHHB, ESTAK, EVDW 
         
         INTEGER :: I, J           !Iteration variables - residue indices
         INTEGER :: K, L           !Atom indices 
         INTEGER :: TYPEI, TYPEJ   !Type of residue (RNA, DNA, protein, ...)
         INTEGER :: TI, TJ         !ID of residue (A, G, ...)
         INTEGER :: TK, TL         !Atom type
         REAL(KIND = REAL64) :: A(3), DA2, DF, DX(3), DCORR
         REAL(KIND = REAL64) :: THIS_EHB, THIS_ESTAK, THIS_EVDW, NB
         LOGICAL :: HBEXIST

         INTEGER :: STACKUNIT, HBUNIT, EXCLVUNIT

         STACKUNIT = GETUNIT()

#ifdef FOR_ANALYSIS
         OPEN(STACKUNIT, FILE="Dbg_Stacking.dat", STATUS='UNKNOWN')
         HBUNIT = GETUNIT()
         OPEN(HBUNIT, FILE="Dbg_Hbonding_E.dat", STATUS='UNKNOWN')
         EXCLVUNIT = GETUNIT()
         OPEN(EXCLVUNIT, FILE="Dbg_ExcludedV.dat", STATUS='UNKNOWN')
#endif    
         EHHB = 0.0D0
         ESTAK = 0.0D0
         EVDW = 0.0D0
         F(1:NOPT) = 0.0D0
         BOCC(1:NRES) = 0

         DO I=1,NRES-1
            DO J=I+1,NRES
               !QUERY: Why are we using first to last here?
               ! Wouldn't it be better to use the distance first to first or last to last?
               !K = BLIST(I)
               !L = BLIST(J-1)+1
               !replaced Blist with RESSTART and RESFINAL to get first and last atom id for all res
               
               K = RESFINAL(I)
               L = RESSTART(J)
               A(1:3) = X(3*K-2:3*K) - X(3*L-2:3*L)
               DA2 = DOT_PRODUCT(A,A)
               !QUERY: This really should be a variable, not a magic number!
               !Check residues are close enough for interactions
               IF (DA2 .GT. NBCUT) THEN
                  CYCLE
               ENDIF
               !Now calculate interactions between residues
               TYPEI = RESTYPES(I)
               TYPEJ = RESTYPES(J)
               TI = BTYPE(I)
               TJ = BTYPE(J)

               !If restype is 0 it is RNA and if it is 1 it is DNA
               !Currently we can do RNA-RNA or DNA-DNA interactions
               IF (((TYPEI.EQ.0).OR.(TYPEI.EQ.1)).AND.((TYPEJ.EQ.0).OR.(TYPEJ.EQ.1))) THEN
                  !Hydrogen bonding between nucleotides
                  CALL ENERGY_HB(I, J, TYPEI, TYPEJ, NOPT, X, F, THIS_EHB, HBEXIST)
                  EHHB = EHHB + THIS_EHB                       
                  !Stacking energy
                  !CALL NA_STACKV(NOPT, BLIST(I),BLIST(J),TI,TJ,F,X,THIS_ESTAK,STACKUNIT)  !old Stacking
                  CALL NA_STACKV2(NOPT, BLIST(I),BLIST(J),TI,TJ,F,X,THIS_ESTAK)
                  ESTAK = ESTAK + THIS_ESTAK
#ifdef FOR_ANALYSIS
                  IF (ABS(THIS_EHB) .GT. 0.3D0) THEN
                     WRITE(HBUNIT,'(4I6,2F15.7)') I, J, TI, TJ, DSQRT(DA2), THIS_EHB
                  ENDIF
#endif
              
               !If restype is 2 for both this is protein-protein
               ELSEIF ((TYPEI.EQ.2).AND.(TYPEJ.EQ.2)) THEN 
                  CYCLE
               ELSE
                  CYCLE
               ENDIF
               !Now calculate excluded volume interactions
               DO K = RESSTART(I),RESFINAL(I)
                  DO L = RESSTART(J),RESFINAL(J)
                     TK = IAC(K)
                     TL = IAC(L)
                     !make sure we skip the bonded particles for neighbouring res
                     IF ((J-I).EQ.1) THEN
                        IF ((TYPEI.EQ.0).AND.(TYPEJ.EQ.0)) THEN
                           !For RNA ignore O3-P
                           IF ((TK.EQ.1).AND.(TL.EQ.3)) CYCLE
                        ELSEIF ((TYPEI.EQ.1).AND.(TYPEJ.EQ.1)) THEN
                           !For DNA do the same
                           IF ((TK.EQ.1).AND.(TL.EQ.3)) CYCLE
                        ENDIF                  
                     ENDIF
                     !additionally skip neighbouring nucleotide bb/sugar to bb/sugar distances (types 1,2,3,4)
                     IF ((J-I).EQ.1) THEN
                        IF ((TYPEI.EQ.0).AND.(TYPEJ.EQ.0)) THEN
                           IF ((TK.LE.7).AND.(TL.LE.7)) CYCLE
                        END IF
                     ENDIF
                     A(1:3) = X(3*K-2:3*K) - X(3*L-2:3*L)
                     DA2 = DOT_PRODUCT(A,A) 
                     !Skip if the distance is too large
                     IF (DA2.GT.RCUT2_EXCLV) CYCLE
                     DCORR = 0.0D0
                     IF (ABS(I-J).EQ.1) DCORR = 0.8D0
                     CALL ENERGY_EXV(NOPT, X, F, K, L, TK, TL, DA2, DCORR, THIS_EVDW)
                      
                     EVDW = EVDW + THIS_EVDW                       
#if FOR_ANALYSIS
                     WRITE(EXCLVUNIT,'(4I6,4F15.7)') I, J , K, L, DSQRT(DA2), NBCT2(TK,TL), THIS_EVDW
#endif
   
                  ENDDO
               ENDDO                
            ENDDO
         ENDDO
#ifdef FOR_ANALYSIS
         CLOSE(HBUNIT)
         CLOSE(STACKUNIT)
         CLOSE(EXCLVUNIT)
#endif  
      END SUBROUTINE E_NONBONDED
   
END MODULE MOD_NONBONDED
