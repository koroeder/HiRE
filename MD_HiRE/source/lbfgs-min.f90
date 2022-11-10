MODULE MINIMISATION
   USE MD_COMMONS, ONLY: NATOMS, MYUNIT, EPOT
   USE NUMKIND
   IMPLICIT NONE
   REAL(KIND = REAL64), ALLOCATABLE, DIMENSION(:) :: DIAG, W
   LOGICAL :: STUCK
   INTEGER :: NQTOT, ITER, POINT
   INTEGER :: ISPT, IYPT, NPT
   ! keyword set inputs for minimiser
   ! history size
   INTEGER :: MUPDATE=4
   ! maximum number of iterations
   INTEGER :: INTMAX=10000
   ! convergence limit
   REAL(KIND = REAL64) :: EPS=1.0D-5
   ! cold fusion limit
   REAL(KIND = REAL64) :: COLDFUSIONLIMIT = -1.0D6
   ! cold fusion hit?
   LOGICAL :: COLDFUSION=.FALSE.
   ! maximum energy increase
   REAL(KIND = REAL64) :: MAXERISE = 0.1D-3
   !Initial guess for Hessian elements
   REAL(KIND = REAL64) :: DGUESS = 0.1D0
   !Maximum decrease in energy allowed
   REAL(KIND = REAL64) :: MAXEFALL = -HUGE(1.0D0)
   REAL(KIND = REAL64) :: RMS = 0.0D0                   !Current RMS force

   CONTAINS

      SUBROUTINE MINIMISE(X,ENERGY,ITDONE,MFLAG)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(INOUT) :: X(3*NATOMS)
         REAL(KIND=REAL64), INTENT(OUT) :: ENERGY
         INTEGER, INTENT(OUT) :: ITDONE
         LOGICAL, INTENT(OUT) :: MFLAG

         CALL MINLBFGS(3*NATOMS, MUPDATE, X, MFLAG, ENERGY, ITDONE, .TRUE.) 
         CALL DEALLOC_LBFGS()
      END SUBROUTINE MINIMISE

      SUBROUTINE ALLOC_LBFGS(N,M)
         INTEGER, INTENT(IN) :: N, M
         CALL DEALLOC_LBFGS()
         ALLOCATE(DIAG(N))       
         ALLOCATE(W(N*(2*M+1)+2*M))         
      END SUBROUTINE ALLOC_LBFGS

      SUBROUTINE DEALLOC_LBFGS()
         IF (ALLOCATED(DIAG)) DEALLOCATE(DIAG)
         IF (ALLOCATED(W)) DEALLOCATE(W)
      END SUBROUTINE DEALLOC_LBFGS

      ! modified L-BFGS algorithm Jorge Nocedal 1990, modifications D J Wales
      SUBROUTINE MINLBFGS(N,M,XCOORDS,MFLAG,ENERGY,ITDONE)
         USE PORFUNCS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD
         IMPLICIT NONE
         INTEGER, INTENT(IN)                  :: N
         INTEGER, INTENT(IN)                  :: M
         INTEGER, INTENT(OUT)                 :: ITDONE
         REAL(KIND = REAL64), INTENT(INOUT)   :: XCOORDS(3*NATOMS)      
         REAL(KIND = REAL64), INTENT(OUT)     :: ENERGY
         LOGICAL, INTENT(INOUT)               :: MFLAG

         REAL(KIND = REAL64) :: GRAD(3*NATOMS), WTEMP(3*NATOMS), XSAVE(N), GNEW(3*NATOMS)
         REAL(KIND = REAL64) :: DUMMY, ENEW, GNORM, STP, YS, YY, SQ, YR, BETA
         REAL(KIND = REAL64) :: DOT1, DOT2, OVERLAP, SLENGTH, DDOT
         INTEGER :: J1, BOUND, CP, INMC, IYCN, ISCN, NFAIL, NDECREASE

         CALL ALLOC_LBFGS(N,M)

         NFAIL=0
         ITER=0
         ITDONE=0

         CALL HIRE_ENERGY_GRAD(3*NATOMS, XCOORDS, ENERGY, GRAD)

         !  Catch cold fusion  and discard.
         IF (ENERGY.LT.COLDFUSIONLIMIT) THEN
            WRITE(MYUNIT,'(A,G20.10)') 'ENERGY=',ENERGY
            WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENERGY,COLDFUSIONLIMIT
            ENERGY=1.0D6
            EPOT=1.0D6
            RMS=1.0D0         
            COLDFUSION=.TRUE.  ! set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
            RETURN
         ENDIF

         EPOT=ENERGY

         IF (MOD(ITDONE,10).EQ.0) THEN
            WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
         END IF

         !  Termination test. 
10       MFLAG=.FALSE.
         IF (RMS.LE.EPS) THEN 
            MFLAG=.TRUE.
               WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
            RETURN
         ENDIF
         
         IF (ITDONE.EQ.ITMAX) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
            RETURN
         ENDIF


         IF (ITER.EQ.0) THEN
            IF (N.LE.0.OR.M.LE.0) THEN
               WRITE(MYUNIT,'(A)') ' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)'
               STOP
            ENDIF
            POINT=0
            MFLAG=.FALSE.
            DO J1=1,N
               DIAG(J1)=DGUESS
            ENDDO

            !     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
            !     ---------------------------------------
            !     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
            !         OTHER TEMPORARY INFORMATION.
            !     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
            !     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
            !         IN THE FORMULA THAT COMPUTES H*G.
            !     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
            !         STEPS.
            !     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
            !         GRADIENT DIFFERENCES.
            !
            !     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
            !     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.

            ISPT= N+2*M    ! index for storage of search steps
            IYPT= ISPT+N*M ! index for storage of gradient differences

            ! NR step for diagonal inverse Hessian
            DO J1=1,N
               DUMMY=-GRAD(J1)*DIAG(J1)
               W(ISPT+J1)=DUMMY
               W(J1)=DUMMY
            ENDDO
            GNORM=DSQRT(DDOT(N,GRAD,1,GRAD,1))
            
            ! Make the first guess for the step length cautious.
            STP=MIN(1.0D0/GNORM,GNORM)

         ELSE 
            BOUND=ITER
            IF (ITER.GT.M) BOUND=M
            YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
            ! Update estimate of diagonal inverse Hessian elements
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (ABS(YY).LT.1.0D-50) THEN
               YY=SIGN(1.0D-50,YY)
            ENDIF
            IF (ABS(YS).LT.1.0D-50) THEN
               YS=SIGN(1.0D-50,YS)
            ENDIF
            DO J1=1,N
               DIAG(J1)= YS/YY
            ENDDO

            !     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
            !     "Updating quasi-Newton matrices with limited storage",
            !     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
            !     ---------------------------------------------------------
            CP= POINT
            IF (POINT.EQ.0) CP=M
            W(N+CP)= 1.0D0/YS
            DO J1=1,N
               W(J1)= -GRAD(J1)
            ENDDO
            CP= POINT
            DO J1= 1,BOUND
               CP=CP-1
               IF (CP.EQ.-1) CP=M-1
               SQ=DDOT(N,W(ISPT+CP*N+1),1,W,1)
               INMC=N+M+CP+1
               IYCN=IYPT+CP*N
               W(INMC)=W(N+CP+1)*SQ
               CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
            ENDDO
         
            DO J1=1,N
               W(J1)=DIAG(J1)*W(J1)
               IF (ABS(W(J1)).GT.1.0D2) W(J1)=SIGN(1.0D2,W(J1)) ! DJW
            ENDDO

            DO J1=1,BOUND
               YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
               BETA= W(N+CP+1)*YR
               INMC=N+M+CP+1
               BETA= W(INMC)-BETA
               ISCN=ISPT+CP*N
               CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
               CP=CP+1
               IF (CP.EQ.M) CP=0
            ENDDO
            STP=1.0D0  
         ENDIF
         !
         !  Store the new search direction
         !
         IF (ITER.GT.0) THEN
            DO J1=1,N
               W(ISPT+POINT*N+J1)= W(J1)
            ENDDO 
         ENDIF
         DOT1=SQRT(DDOT(N,GRAD,1,GRAD,1))

         !  Overflow has occasionally occurred here.
         !  We only need the sign of the overlap, so use a temporary array with
         !  reduced elements.
         DUMMY=1.0D0
         DO J1=1,N
            IF (ABS(W(J1)).GT.DUMMY) DUMMY=ABS(W(J1))
         ENDDO
         DO J1=1,N
            WTEMP(J1)=W(J1)/DUMMY
         ENDDO
         DOT2=SQRT(DDOT(N,WTEMP,1,WTEMP,1))
         OVERLAP=0.0D0
         IF (DOT1*DOT2.NE.0.0D0) THEN
            OVERLAP=DDOT(N,GRAD,1,WTEMP,1)/(DOT1*DOT2)
         ENDIF
         IF (OVERLAP.GT.0.0D0) THEN
            IF (DEBUG) WRITE(MYUNIT,'(A)') 'Search direction has positive projection onto gradient - reversing step'
            DO J1=1,N
               W(ISPT+POINT*N+J1)= -W(J1)  !!! DJW, reverses step
               WTEMP(J1)= -W(J1)/DUMMY  !!! DJW, reverses step
            ENDDO
         ENDIF

         DO J1=1,N
            W(J1)=GRAD(J1)
         ENDDO
         SLENGTH=0.0D0
         DO J1=1,N
            SLENGTH=SLENGTH+WTEMP(J1)**2
         ENDDO
      
         ! Taking out the magnitude of the largest element can prevent overflow.
         SLENGTH=SQRT(SLENGTH)*DUMMY
         IF (STP*SLENGTH.GT.MAXBFGS) THEN
            STP=MAXBFGS/SLENGTH
         ENDIF

         ! We now have the proposed step.
         ! Save XCOORDS here so that we can undo the step reliably 
         XSAVE(1:N)=XCOORDS(1:N)
         DO J1=1,N
            XCOORDS(J1)=XCOORDS(J1)+STP*W(ISPT+POINT*N+J1)
         ENDDO

         NDECREASE=0

20       CALL HIRE_ENERGY_GRAD(3*NATOMS, XCOORDS, ENERGY, GRAD)

         !  Catch cold fusion and discard.
         IF (ENEW.LT.COLDFUSIONLIMIT) THEN
            WRITE(MYUNIT,'(A,2G20.10)') ' Cold fusion diagnosed - step discarded; energy and threshold=',ENEW,COLDFUSIONLIMIT
            ENERGY=1.0D6
            ENEW=1.0D6
            POTEL=1.0D6
            RMS=1.0D0
            COLDFUSION=.TRUE.  ! set COLDFUSION=.TRUE. so that ATEST=.FALSE. in MC
            RETURN
         ENDIF

         IF (((ENEW-ENERGY.LE.MAXERISE)).AND.(ENEW-ENERGY.GT.MAXEFALL)) THEN
            ITER=ITER+1
            ITDONE=ITDONE+1
            ENERGY=ENEW
            DO J1=1,3*NATOMS
               GRAD(J1)=GNEW(J1)
            ENDDO
            IF (DEBUG) THEN
               WRITE(MYUNIT,'(A,G20.10,G20.10,A,I6,A,F13.10)') ' Energy and RMS force=',ENERGY,RMS,' after ',ITDONE, &
                                                               ' LBFGS steps, step:',STP*SLENGTH
            ENDIF

         !  May want to prevent the PE from falling too much if we are trying to visit all the
         !  PE bins. Halve the step size until the energy change is in range.
         ELSEIF (ENEW-ENERGY.LE.MAXEFALL) THEN
            IF (NDECREASE.GT.15) THEN
               NFAIL=NFAIL+1
               WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find an energy in the required range, NFAIL=',NFAIL           
               ! Resetting to XSAVE should be the same as subtracting the step.
               IF(DEBUG) write(*,*) "mymylbfgs> Resetting to saved coords after failed step"
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
               ITER=0   !  try resetting
               IF (NFAIL.GT.20) THEN
                  WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
                  RETURN
               ENDIF
               GOTO 30
            ENDIF
            ! Resetting to XSAVE and adding half the step should be the same as subtracting 
            ! half the step.
            IF(DEBUG) WRITE(*,*) "mymylbfgs> Resetting to saved coords after failed step"
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.5*STP*W(ISPT+POINT*N+J1)
            ENDDO
            STP=STP/2.0D0
            NDECREASE=NDECREASE+1
            IF (DEBUG) WRITE(MYUNIT,'(A,F19.10,A,F16.10,A,F15.8)') &
                            ' energy increased too much from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
            GOTO 20
         ELSE
            ! Energy increased - try again with a smaller step size
            IF (NDECREASE.GT.10) THEN ! DJW
               NFAIL=NFAIL+1
               WRITE(MYUNIT,'(A,G20.10)') ' in mylbfgs LBFGS step cannot find a lower energy, NFAIL=',NFAIL
               ! Resetting to XSAVE should be the same as subtracting the step. 
               XCOORDS(1:N)=XSAVE(1:N)
               GRAD(1:N)=GNEW(1:N) ! GRAD contains the gradient at the lowest energy point
               ITER=0   !  try resetting
               IF (NFAIL.GT.5) THEN         
                  WRITE(MYUNIT,'(A)') ' Too many failures - giving up '
                  RETURN
               ENDIF
               GOTO 30
            ENDIF
            ! Resetting to XSAVE and adding 0.1 of the step should be the same as subtracting 
            ! 0.9 of the step. 
            XCOORDS(1:N)=XSAVE(1:N)
            DO J1=1,N
               XCOORDS(J1)=XCOORDS(J1)+0.1D0*STP*W(ISPT+POINT*N+J1)
            ENDDO
            STP=STP/1.0D1
            NDECREASE=NDECREASE+1
            IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,G20.10,A,G20.10)') &
                            ' energy increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
            GOTO 20
         ENDIF

   !     Compute the new step and gradient change
   30    NPT=POINT*N
         DO J1=1,N
            W(ISPT+NPT+J1)= STP*W(ISPT+NPT+J1) ! save the step taken
            W(IYPT+NPT+J1)= GRAD(J1)-W(J1)     ! save gradient difference: W(1:N) contains the old gradient
            !  Prevent overflow of dot product for large gradients.
            IF (ABS(W(IYPT+NPT+J1)).GT.1.0D10) W(IYPT+NPT+J1)=SIGN(1.0D10,W(IYPT+NPT+J1))
         ENDDO

         POINT=POINT+1
         IF (POINT.EQ.M) POINT=0
         IF (CENT) CALL CENTRE2(XCOORDS)

         GOTO 10
         RETURN
      END SUBROUTINE MINLBFGS

END MODULE MINIMISATION
