      SUBROUTINE QUENCH(QTEST,NP,ITER,TIME,QDONE,P)
      USE MODHESS  
      USE GENRIGID
#ifdef __SPARSE
      USE MODSPARSEHESS
      USE SHIFT_HESS
      USE INERTIA_MOD
#endif /* __SPARSE */
      USE COMMONS
      USE PREC
      use porfuncs
      IMPLICIT NONE

      INTEGER, INTENT(IN)              :: NP, ITER
      INTEGER, INTENT(OUT)             :: QDONE
      REAL(KIND = REAL64), INTENT(OUT) :: P(3*NATOMS), TIME
      LOGICAL, INTENT(IN)              :: QTEST

      REAL(KIND = REAL64)   :: GRAD(3*NATOMS), EREAL
      LOGICAL               :: CFLAG
      INTEGER               :: NOPT

      INTEGER               :: I, J1, HORDER, NUM_ZERO_EVS

      ! Timing variables
      REAL(KIND = REAL64)   :: T_LBFGS_START, T_LBFGS_END, T_DIAG_START, T_DIAG_END

      ! Sparse/Hessian
      REAL(KIND = REAL64)   :: FRQS(3*NATOMSALLOC), EVALUES(3*NATOMSALLOC)
      REAL(KIND = REAL64)   :: LOG_PROD_EV, LARGEST_ZERO, SMALLEST_NONZERO
      INTEGER               :: ATTEMPTS
      LOGICAL               :: TS_FOUND

      ! FEBH calculation
      REAL(KIND = REAL64)   :: NEWPFMIN, ITDET, MI1, MI2, MI3

      ! the following required to call the LAPACK routine DSYEV
      INTEGER               :: INFO
      INTEGER, PARAMETER    :: LWORK=10000 ! the dimension is set arbitrarily
      REAL(KIND = REAL64)   :: WORK(LWORK)

!  QTEST is set for the final quenches with tighter convergence criteria.
      IF (QTEST) THEN
         GMAX=CQMAX
      ELSE
         GMAX=BQMAX
      ENDIF

      NOPT=3*NATOMS
      QDONE=0
      DO I=1,3*NATOMS
         P(I)=COORDS(I,NP)
      ENDDO

      ! Added timing for LBFGS call
      CALL CPU_TIME(T_LBFGS_START)
      CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.)
      CALL CPU_TIME(T_LBFGS_END)
      IF (DEBUG) WRITE(MYUNIT, '(A, F10.3)') 'quench> Time to minimise (s):', T_LBFGS_END - T_LBFGS_START

      ! Set parameters that may not be initialised if we aren't doing FEBH etc.
      FEBH_POT_ENE=EREAL
      HORDER=-1
      LOG_PROD_EV=-1.0D0

!###############################LEFT OFF HERE


      IF (FEBHT .AND. CFLAG) THEN
         ! Calculate the free energy
         NUM_ZERO_EVS=6
         IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
         ALLOCATE(HESS(3*NATOMS,3*NATOMS))
         CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
         CALL MASSWT()
         
#ifdef __SPARSE
         IF (SPARSET) THEN
             IF (DEBUG) THEN
                WRITE(MYUNIT, *) 'quench> Using sparse cholesky decomposition to calculate determinant.'
             END IF
             SHIFT_ARRAY(:) = 1.0D0
             CALL CPU_TIME(T_DIAG_START)
             CALL SHIFT_HESS_ZEROS(P, SHIFT_ARRAY)
             CALL FILTER_ZEROS(HESS, ZERO_THRESH)
             CALL GET_DETERMINANT(3*NATOMS, LOG_PROD_EV, NQ(NP))
             CALL CPU_TIME(T_DIAG_END)
             IF (DEBUG) WRITE(MYUNIT,'(A, F10.3)') 'quench> Time to diagonalise with sparse routines (s):',T_DIAG_END-T_DIAG_START
             IF (LOG_PROD_EV < SMALL_DOUBLE) THEN
                WRITE(MYUNIT, '(A)') 'quench> Log product from SuiteSparse is -infinity. The hessian ' &
     &                               'is probably not positive definite. If using a threshold ' &
     &                               'on zero values, try reducing it.'

                ! Changed to convergence failure so that the run can continue, as for ts location. DJW
                CFLAG=.FALSE.
                RETURN
             END IF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For debugging only, this *should* be very slow for cases where we use the sparse approach: 
!             PRINT *, "Sparse:", LOG_PROD_EV
!             CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
!             CALL MASSWT()
!             CALL SHIFT_HESS_ZEROS(P, (/ 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0, 1.0D0 /))
!             CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)             
!             PRINT *, "DSYEV:", SUM(LOG(EVALUES))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ELSE
#endif /* __SPARSE */
            CALL CPU_TIME(T_DIAG_START)
            CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)
            CALL CPU_TIME(T_DIAG_END)
            IF (DEBUG) THEN
                WRITE(MYUNIT, '(A, F10.3)') 'quench> Time to diagonalise with DSYEV (s):', T_DIAG_END - T_DIAG_START
                WRITE(MYUNIT, '(A)') "quench> Eigenvalues:"
                NUM_ZERO_EVS=6
                WRITE(MYUNIT, '(3F12.3)') EVALUES(1:3*NATOMS)
                WRITE(MYUNIT, '(A)') "quench> Frequencies in wavenumbers:"
                DO J1=1,3*NATOMS
                  FRQS(J1)=108.575D0*SQRT(ABS(EVALUES(J1)))
                ENDDO
                WRITE(MYUNIT, '(3F12.3)') FRQS(1:3*NATOMS)
            END IF
! Check that we don't have any proper negative eigenvalues (the magnitude should drop).
! We do it this way, since some of the zeros come out negative as well, so you can't 
! just check the sign.
            TS_FOUND = (ABS(EVALUES(1)) > ABS(EVALUES(NUM_ZERO_EVS + 1)))
            IF (TS_FOUND) THEN
                WRITE(MYUNIT, '(A, I10, A)') 'quench> Quench ', NQ(NP), ' converged to a transition state.'
                IF (DEBUG) THEN
                   WRITE(MYUNIT, '(A)') 'quench> Eigenvalues'
                   WRITE(MYUNIT, '(A)') '======================='
                   WRITE(MYUNIT, '(6F20.12)') EVALUES(1:(3+NUM_ZERO_EVS))
                   WRITE(MYUNIT, '(A)') '======================='
                END IF
                CFLAG=.FALSE.
                RETURN                           
            END IF
            LOG_PROD_EV=SUM(DLOG(EVALUES(NUM_ZERO_EVS + 1:3*NATOMS)))
            IF ((MIN_ZERO_SEP .GT. 0.0D0) .AND. (MAX_ATTEMPTS .GT. 0)) THEN
                ! Square the minimum zero separation, since this refers to separation of frequencies.
                DO ATTEMPTS = 1, MAX_ATTEMPTS            
                    ! Identify the largest zero and smallest non-zero
                    EVALUES(:) = ABS(EVALUES(:))
                    LARGEST_ZERO = 0.D0
                    LARGEST_ZERO = MAX(LARGEST_ZERO, MAXVAL(EVALUES(1:NUM_ZERO_EVS)))
                    SMALLEST_NONZERO = 1.0D100
                    SMALLEST_NONZERO = MIN(SMALLEST_NONZERO, MINVAL(EVALUES(NUM_ZERO_EVS+1:)))
                    ! If the separation of zeros and non-zeros is too small, reduce the convergence
                    ! threshold by an order of magnitude, along with the corresponding sloppy or tight
                    ! convergence threshold (for future quenches).
                    IF ((SMALLEST_NONZERO / LARGEST_ZERO) < MIN_ZERO_SEP) THEN
                        WRITE(MYUNIT, '(A,I8)') 'quench> Attempt: ', ATTEMPTS
                        WRITE(MYUNIT, '(A,E15.7)') 'quench> Current separation is ', (SMALLEST_NONZERO / LARGEST_ZERO)
                        WRITE(MYUNIT, '(A,E15.7)') 'quench> Target separation is ', MIN_ZERO_SEP
                        GMAX = GMAX * 1.0D-1
                        IF (QTEST) THEN
                            CQMAX = GMAX
                            WRITE(MYUNIT, '(A,F12.8)') 'quench> Lowering tight convergence to ', CQMAX
                        ELSE
                            BQMAX = GMAX
                            WRITE(MYUNIT, '(A,F12.8)') 'quench> Lowering sloppy convergence to ', BQMAX
                        END IF
                        CALL MYLBFGS(NOPT,MUPDATE,P,.FALSE.,GMAX,CFLAG,EREAL,MAXIT,ITER,.TRUE.)
                        CALL POTENTIAL(P,GRAD,EREAL,.TRUE.,.TRUE.)
                        CALL MASSWT()
                        CALL DSYEV('N','L',3*NATOMS,HESS,3*NATOMS,EVALUES,WORK,LWORK,INFO)
                        EVALUES(:) = ABS(EVALUES(:))
                        ! Recalculate separation.
                        LARGEST_ZERO = 0.D0
                        LARGEST_ZERO = MAX(LARGEST_ZERO, MAXVAL(EVALUES(1:NUM_ZERO_EVS)))
                        SMALLEST_NONZERO = 1.0D100
                        SMALLEST_NONZERO = MIN(SMALLEST_NONZERO, MINVAL(EVALUES(NUM_ZERO_EVS+1:)))
                    END IF
                    IF ((SMALLEST_NONZERO / LARGEST_ZERO) >= MIN_ZERO_SEP) THEN
                        ! If we have converged properly, then we're ok.
                        IF (DEBUG) THEN
                           WRITE(MYUNIT, '(A,E15.7)') 'quench> Converged. Separation of zeros: ', (SMALLEST_NONZERO / LARGEST_ZERO)
                        END IF
                        EXIT
                    END IF
                    IF (ATTEMPTS == MAX_ATTEMPTS) THEN
                        WRITE(MYUNIT, '(A)') 'quench> Failed to achieve desired separation of zeros and non-zeros.'
                        WRITE(MYUNIT, '(A)') 'quench> Magnitudes of lowest eigenvalues:'
                        WRITE(MYUNIT, '(F20.12)') EVALUES(1:NUM_ZERO_EVS + 6)
                        STOP 'Cannot achieve desired separation of zeros and non-zeros. Please check your input.'
                    END IF
                END DO
                LOG_PROD_EV = SUM(DLOG(EVALUES(NUM_ZERO_EVS+1:)))
            END IF
#ifdef __SPARSE
         END IF
#endif /* __SPARSE */

! Continue with FEBH calculation
         ITDET=0.0D0
         !Assume we have C1 symmetry in this case
         HORDER=1
         IF (USEROT) THEN
            CALL MYINERTIA(P,NATOMS,ITDET,MI1,MI2,MI3)
            PRINT '(A,3G20.10)','MI1,MI2,MI3=',MI1,MI2,MI3
            ITDET=0.0D0
            MI1=-1.0D0; MI2=-1.0D0; MI3=-1.0D0
         ENDIF
         IF (DEBUG) THEN
            IF (QTEST) WRITE(MYUNIT,'(A,I4)') 'quench> Order of the point group of minimum=', HORDER
         ENDIF
         IF (ITDET.NE.0.0D0) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'quench> Inertia determinant=',ITDET
            ITDET=LOG(ITDET)/2.0D0 
         ENDIF

         FEBH_POT_ENE=EREAL ! saves potential energy.

         ! NEWPFMIN ln( n! qrot * qvib * exp(-V/kT) / o ) with o the order of the point group
         CALL MYINERTIA(P,NATOMS,ITDET,MI1,MI2,MI3)
         CALL GETPFMIN(FETEMP,NEWPFMIN,NATOMS,MI1,MI2,MI3,LOG_PROD_EV,EREAL,HORDER,EVALUES)
         IF (FETEMP.EQ.0.D0) THEN
            IF (DEBUG) WRITE(MYUNIT, '(A,F20.12)') 'quench> Potential energy=', FEBH_POT_ENE
            EREAL=EREAL  ! molecule monomers may need the HORDERMIN term and rotation/vibration
         ELSE
            EREAL=-FETEMP*NEWPFMIN          
         ENDIF

      ELSE IF (FEBHT .AND. (.NOT. CFLAG)) THEN
         WRITE(MYUNIT, '(A)') 'quench> Quench did not converge, not calculating free energy and adding 1E10 to energy.'
         EREAL=EREAL + 1.0D10
      ENDIF

      POTEL=EREAL

      IF (CFLAG) QDONE=1
      IF (.NOT.CFLAG) THEN
         IF (QTEST) THEN
            WRITE(MYUNIT,'(A,I6,A)') 'quench> WARNING - Final Quench ',NQ(NP),'  did not converge'
         ELSE
            WRITE(MYUNIT,'(A,I7,A)') 'quench> WARNING - Quench ',NQ(NP),'  did not converge'            
         ENDIF
      ENDIF

      CALL MYCPU_TIME(TIME)

      ! Write free energies to an output file after the quench.
      ! Quench=NQ(NP) 
      ! Potential energy=FEBH_POT_ENE
      ! Free energy=EREAL
      ! Markov energy=EPREV(NP)
      ! Harmonic superposition contribution=EREAL - FEBH_POT_ENE
      ! Time=TIME
      IF (FEBHT.AND.(MOD(NQ(NP)-1,PRTFRQ).EQ.0)) THEN
         IF (.NOT. QTEST) THEN
            WRITE(FE_FILE_UNIT, '(I12,8X,4(G20.12,2X),F18.1,2X)')         NQ(NP), FEBH_POT_ENE, EREAL - FEBH_POT_ENE,&
     &                                                               EREAL, EPREV(NP), TIME
         ELSE
            WRITE(FE_FILE_UNIT, '(A1,I11,8X,4(G20.12,2X),F18.1,2X)') 'F', NQ(NP), FEBH_POT_ENE, EREAL - FEBH_POT_ENE,&
     &                                                               EREAL, EPREV(NP), TIME
         END IF
      END IF


      ! SAVEQ is .TRUE. for quenches. For final quenches, it 
      ! is set to .FALSE. (in finalq.f) and so the checks are skipped.
      IF (SAVEQ) THEN        
         CALL GSAVEIT(EREAL,P,NP)
         IF (NSAVE==0) QMIN(1)=min(QMIN(1),EREAL)
      ENDIF

      DO J1=1,3*NATOMS
         COORDS(J1,NP)=P(J1)
      ENDDO

      RETURN
      END SUBROUTINE QUENCH
