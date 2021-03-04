!  Final quenching of the lowest saved minima to the tight convergence criterion specified

      SUBROUTINE FINALQ
      USE COMMONS
      USE PREC
      USE GENRIGID

      IMPLICIT NONE

      INTEGER :: J1, J2, NP, ITERATIONS, QDONE
      REAL(KIND = REAL64) :: TIME, SCREENC(3*NATOMSALLOC)

      SAVEQ=.FALSE.
      NQ(:)=0
      MAXIT=MAXIT2
      IF (MPIT) THEN
         NP=MYNODE+1
      ELSE
         NP=1
      ENDIF

      ! Adding in clarification of what is being tightly quenched
      WRITE(MYUNIT,'(A)') 'Tightly converging the SAVE lowest energy minima found'
      WRITE(MYUNIT,'(A)') 'NOTE: these may NOT match the other output files - see below for a sorted list of Lowest minima'
      DO J1=1,NSAVE
         IF (QMIN(J1).LT.1.0D10) THEN
            DO J2=1,3*NATOMS
               COORDS(J2,NP)=QMINP(J1,J2)
            ENDDO
            NQ(NP)=NQ(NP)+1
            IF (RELAXFQ.OR.HYBRIDMINT) THEN
               RIGIDINIT = .FALSE.
            ENDIF
            IF(UNFREEZEFINALQ) FROZEN(:)=.FALSE.  ! unfreeze everything before the final quenches
            CALL QUENCH(.TRUE.,NP,ITERATIONS,TIME,QDONE,SCREENC)
            WRITE(MYUNIT,'(A,I6,A,G20.10,A,I5,A,G15.7,A,F12.2)') 'Final Quench ',NQ(NP),' energy=', &
     &                POTEL,' steps=',ITERATIONS,' RMS force=',RMS,' time=',TIME-TSTART
            QMIN(J1)=POTEL
            DO J2=1,3*NATOMS
               QMINP(J1,J2)=COORDS(J2,NP)
            ENDDO
         ENDIF
      ENDDO

      NSAVE=MIN(NSAVE,NQ(NP))
      ! Re-sort the saved minima now that they have been tightly converged
      CALL GSORT2()
      ! Print a list of sorted minima energies which will 100% match other output files
      WRITE(MYUNIT,'(A)') 'After re-sorting, the lowest found minima are (lowest free energy subtracted if applicable):'
      DO J1=1,NSAVE
         WRITE(MYUNIT,'(A,I6,A,G20.10)') 'Lowest Minimum ',J1,' Energy= ',QMIN(J1)
      ENDDO

      RETURN
      END
