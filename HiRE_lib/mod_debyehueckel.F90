!> @file
!> Contains MOD_DEBYEHUECKEL to deal with Debye-Hueckel contribution

!> Module contains all routines and variables to calculate the Debye-Hueckel contribution 
MODULE MOD_DEBYEHUECKEL
   USE PREC_HIRE
   USE NBDEFS
   IMPLICIT NONE
   !> Charged particles ids
   INTEGER, ALLOCATABLE :: CHARGED_IDS(:)
   !> NUmber of charged particles
   INTEGER :: NCHARGED = 0
   CONTAINS
      !> Routine to allocate all required arrays
      SUBROUTINE INIT_DH()
         USE NAPARAMS, ONLY: SCORE_RNA, DL
         DL = SCORE_RNA(16)
         CALL FIND_CHARGED_PARTICLES()
      END SUBROUTINE INIT_DH

      !> Routine to set array of ids for charged particles (the routine can be called at any time)
      SUBROUTINE FIND_CHARGED_PARTICLES()
         USE VAR_DEFS, ONLY: CHATM, NPARTICLES
         IMPLICIT NONE
         INTEGER :: IDX
         INTEGER :: I
         REAL(KIND = REAL64), PARAMETER :: EPSZERO = 1.0D-6

         CALL DEALLOC_DH()
         ! first find all charged particles
         NCHARGED = 0
         DO I=1,NPARTICLES
            IF (ABS(CHATM(I)).GT.EPSZERO) THEN
               NCHARGED = NCHARGED + 1
            END IF
         END DO
         ALLOCATE(CHARGED_IDS(NCHARGED))
         ! now add the indices to the array
         IDX = 0
         DO I=1,NPARTICLES
            IF (ABS(CHATM(I)).GT.EPSZERO) THEN
               IDX = IDX + 1
               CHARGED_IDS(IDX) = I
            END IF
         END DO
         !WRITE(*,*) " find_charged> Found ", NCHARGED, " particles"
      END SUBROUTINE FIND_CHARGED_PARTICLES

      SUBROUTINE DEALLOC_DH()
         IMPLICIT NONE
         IF (ALLOCATED(CHARGED_IDS)) DEALLOCATE(CHARGED_IDS)
      END SUBROUTINE DEALLOC_DH  

      !> Subroutine only iterating over charged particles, increasing computational efficiency
      SUBROUTINE DH_ENERGY(NOPT, X, F, EDH)
         USE VAR_DEFS, ONLY: CHATM
         USE NAPARAMS, ONLY: DHCUT
#if FOR_ANALYSIS        
         USE UTILS_IO, ONLY: GETUNIT
#endif         
         IMPLICIT NONE
     
         INTEGER, INTENT(IN) :: NOPT                   !should be 3*NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(NOPT)    !input coordinates
         REAL(KIND = REAL64), INTENT(OUT) :: F(NOPT)   !force from bonds
         REAL(KIND = REAL64), INTENT(OUT) :: EDH

         INTEGER :: I, J, IDX, JDX
         REAL(KIND = REAL64) :: RIJ(3), R, CHRGI, CHRGJ, EDHPAIR, DFPAIR, DX(3)

#if FOR_ANALYSIS
         INTEGER :: DHUNIT

         DHUNIT = GETUNIT()
         OPEN(DHUNIT, FILE="Dbg_DebyeHueckel.dat", STATUS='UNKNOWN')
#endif

         EDH = 0.0D0
         F(1:NOPT) = 0.0D0
         
         ! iterate over the charged particle list
         DO IDX=1,NCHARGED-1
            DO JDX=IDX+1,NCHARGED
               I = CHARGED_IDS(IDX)
               J = CHARGED_IDS(JDX)
               ! get charges
               CHRGI = CHATM(I)
               CHRGJ = CHATM(J)
               ! calculate the distance between the atoms
               RIJ(1:3) = X(I*3-2:I*3) - X(3*J-2:3*J)
               R = DSQRT(DOT_PRODUCT(RIJ,RIJ))
               ! use cutoff here, at 35 Angstrom the contributions are smaller than 10^-6
               ! using cutoffs of 50 Angstorm and larger seem to work well cumulatively
               IF (R.LT.DHCUT) THEN
                  ! get the energy for this pair
                  EDHPAIR = 0.0D0
                  DFPAIR = 0.0D0
                  CALL DH_PAIR(R, EDHPAIR, DFPAIR, CHRGI, CHRGJ) 
                  EDH = EDH + EDHPAIR              
                  DX(1:3) = DFPAIR * RIJ(1:3)
                  !WRITE(*,'(2I8,4F15.7)') I, J, R, EDHPAIR, NB, EDHPAIR*NB
#ifdef FOR_ANALYSIS
                  WRITE(DHUNIT,'(2I6,2F7.3,4F15.7)') I, J, CHRGI, CHRGJ, R, EDHPAIR
#endif
                  F((I*3-2):I*3) = F((I*3-2):I*3) - DX(1:3)
                  F((J*3-2):J*3) = F((J*3-2):J*3) + DX(1:3)    
               END IF           
            END DO
         END DO
#ifdef FOR_ANALYSIS
         CLOSE(DHUNIT)
#endif
      END SUBROUTINE DH_ENERGY

      !> DH contribution for a pair of CG particles
      SUBROUTINE DH_PAIR(R, EDH, DF, QI, QJ)
         USE NAPARAMS, ONLY: DIEL, DL
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

