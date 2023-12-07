
MODULE MCmod
  
  USE COMMONS
  USE DEFS_MCRUNS
  USE PREC
  IMPLICIT NONE

  REAL(KIND = REAL64), ALLOCATABLE :: BESTCOORDS(:,:) , SAVECOORDS(:), TEMPCOORDS(:), SCREENC(:)

  CONTAINS

      SUBROUTINE MC(NSTEPS,SCALEFAC)
      USE TWIST_MOD, ONLY: TWISTT
      USE GENRIGID
      USE MOVES
      USE GROUPROTMOD
      USE BP_MOVES_MOD
      USE HIRE_INTERFACE, ONLY: HIRE_SAXS_FORCE
      USE STOCH_FORCE_STEPS, ONLY: STOCHFORCET, GRADMOD_STEP, DOGRADMODSTEP
      USE porfuncs


      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NSTEPS
      REAL(KIND = REAL64), INTENT(IN) :: SCALEFAC
#ifdef MPI
      INCLUDE 'mpif.h'
      INTEGER MPIERR, ITRAJ, NEACCEPT, NTOT
      LOGICAL MPI_FINT
      REAL(KIND = REAL64) :: CTE, TEMPTRAJ(0:NPAR-1), T, BETA(0:NPAR-1)
#endif
      REAL(KIND = REAL64), PARAMETER :: PI=3.141592654D0

      INTEGER             :: J1, J2, J3, JP, REJSTREAK, NDONE, ITERATIONS, QDONE
      INTEGER             :: NSUCCESS(NPAR), NFAIL(NPAR), NFAILT(NPAR), NSUCCESST(NPAR), JBEST(NPAR)
      INTEGER             :: NULLMOVES(NPAR), NULLMOVEST(NPAR)

      REAL(KIND = REAL64) :: GRAD(3*NATOMSALLOC)
      REAL(KIND = REAL64) :: EBEST(NPAR)
      REAL(KIND = REAL64) :: TIME, RANDOM, OPOTEL, MCTEMP 
      INTEGER             :: RNDSEED

      LOGICAL             :: LOPEN, COMPLETE, ATEST, MOVEDT

      CHARACTER(len=10)   :: DATECHAR,TIMECHAR,ZONECHAR
      INTEGER             :: VALUES(8)

      LOGICAL             :: LOOSEFT, LOOSETT
      INTEGER             :: NUCF, NUCT

      REAL(KIND = REAL64) :: SAXSSTEPSIZE, TSAXS1, TSAXS2, ESAXS, SAXSFORCE(3*NATOMS)

      IF (.NOT.ALLOCATED(BESTCOORDS)) ALLOCATE(BESTCOORDS(3*NATOMSALLOC,NPAR))
      IF (.NOT.ALLOCATED(SAVECOORDS)) ALLOCATE(SAVECOORDS(3*NATOMSALLOC))
      IF (.NOT.ALLOCATED(TEMPCOORDS)) ALLOCATE(TEMPCOORDS(3*NATOMSALLOC))
      IF (.NOT.ALLOCATED(SCREENC)) ALLOCATE(SCREENC(3*NATOMSALLOC))
      IF (.NOT.ALLOCATED(TMOVE)) ALLOCATE(TMOVE(NPAR))
      IF (.NOT.ALLOCATED(OMOVE)) ALLOCATE(OMOVE(NPAR))  
      ! initialise RANDOM to 0.0D0 to avoid printing uninitialised variable 
      RANDOM=0.0D0

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'mc> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF

#ifdef MPI
      IF (MPIT) THEN 
         IF (DEBUG) WRITE(MYUNIT,'(A,I6)') 'MPIERR=',MPIERR
         CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
         IF (DEBUG) WRITE(MYUNIT, '(A,2I6)') 'NPAR,MPIERR=',NPAR,MPIERR
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
         JP=MYNODE+1
         ITRAJ=MYNODE
         IF (DEBUG) WRITE(MYUNIT, '(A,3I6)') 'In mc after MPI_MPI_COMM_RANK MPIERR,MYNODE,JP=',MPIERR,MYNODE,JP
      ENDIF
      NEACCEPT=0

#endif

      NDONE=0
      IF (RESTORET) THEN
#ifdef MPI
         CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
#else
         DO JP=1,NPAR
            CALL RESTORESTATE(NDONE,EBEST,BESTCOORDS,JBEST,JP)
         ENDDO
#endif
         WRITE(MYUNIT, '(A,I10)') 'MC> restore NDONE=',NDONE
      ENDIF
      NQ(:)=NDONE

#ifdef MPI
      WRITE(MYUNIT, '(A,I10,A,I10,A)') "Processor", MYNODE+1, " of", NPAR, " speaking:"
      WRITE(MYUNIT, '(A,I10)') 'Number of atoms', NATOMS
#endif

      IF (NACCEPT.EQ.0) NACCEPT=NSTEPS+1
      NQTOT=0

#ifdef MPI
#else
      DO JP=1,NPAR 
#endif
         TMOVE(JP)=.TRUE.
         OMOVE(JP)=.TRUE.
         NSUCCESS(JP)=0
         NFAIL(JP)=0
         NSUCCESST(JP)=0
         NFAILT(JP)=0
#ifdef MPI
#else
      ENDDO
#endif
   
      IF (.NOT.RESTORET) THEN
         IF (SETCENT) CALL SETCENTRE(COORDS)
      ENDIF

!  Calculate the initial energy and save in EPREV
#ifdef MPI
      WRITE(MYUNIT, '(A)')  'Calculating initial energy'

      CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,QDONE,SCREENC)
      NQTOT=NQTOT+1
      WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=', &
     &             POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
      WRITE(MYUNIT, *) 'MC POSTQUENCH'
      IF(DEBUG) WRITE(MYUNIT, '(3F20.10)') COORDS
      CALL FLUSH(MYUNIT)
 
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
      EPREV(JP)=POTEL
      IF (.NOT.RESTORET) EBEST(JP)=POTEL
      BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
      JBEST(JP)=0
      COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)

! Initialisation 

      IF (PTTMIN < 1.0D-6 ) PTTMIN = 1.0D-6 ! to avoid devision by zero
      CTE=(LOG(PTTMAX/PTTMIN))/(NPAR-1)
      CTE=EXP(CTE)

      DO I=0, NPAR-1
         TEMPTRAJ(I)=PTTMIN*CTE**I
         T=TEMPTRAJ(I)
         BETA(I)=1.0D0/T
      ENDDO
      DO I=1, NPAR
         TEMP(I)=TEMPTRAJ(I-1)
      ENDDO
      CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE., .FALSE.)

      WRITE(MYUNIT, '(A, 2G20.10)') 'Temperature range', TEMPTRAJ(0), TEMPTRAJ(NPAR-1)
      WRITE(MYUNIT, '(A, G20.10)') 'For this replica T=', TEMPTRAJ(MYNODE)
      WRITE(MYUNIT, '(A, G20.10)') 'Starting potential energy=', POTEL

      IF (.NOT.RANSEEDT) THEN
         CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
         RNDSEED=values(6)*60 + values(7)
         CALL SDPRND(RNDSEED+MYNODE)
         CALL SDPRND_UNIVERSAL(RNDSEED+NPAR)
         WRITE(MYUNIT, '(A)') ' mc> seeded random number generator with system time'
      ENDIF
      TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
      
#else
      WRITE(MYUNIT,'(A)') 'Calculating initial energy'
      DO JP=1,NPAR
         CALL POTENTIAL(COORDS(:,MYNODE+1),GRAD, POTEL, .TRUE., .FALSE., .FALSE.)          
         WRITE(MYUNIT,'(A,I10)') 'mc calling initial quench'
         CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,QDONE,SCREENC)
         NQTOT=NQTOT+1
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP,']Qu ',NQ(JP),' E=', &
     &             POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ELSE
            WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=', &
     &             POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',POTEL,' t=',TIME-TSTART
         ENDIF
!  EPREV saves the previous energy in the Markov chain.
!  EBEST and JBEST record the lowest energy since the last reseeding and the
!  step it was attained at. BESTCOORDS contains the corresponding coordinates.
 
         EPREV(JP)=POTEL
         IF (.NOT.RESTORET) EBEST(JP)=POTEL
         BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         JBEST(JP)=0
         COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         IF (.NOT.RESTORET) THEN
            EBEST(JP)=POTEL
            BESTCOORDS(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
         ENDIF
      ENDDO
      IF (.NOT.RANSEEDT) THEN
         CALL DATE_AND_TIME(datechar,timechar,zonechar,values)
         RNDSEED=values(6)*60 + values(7)
         CALL SDPRND(RNDSEED)
         CALL SDPRND_UNIVERSAL(RNDSEED+NPAR)
         WRITE(MYUNIT, '(A)') ' mc> seeded random number generator with system time'
      ENDIF
#endif

      IF (NPAR.EQ.1) THEN
         WRITE(MYUNIT,'(A,I10,A)') 'Starting run of ',NSTEPS,' BH steps'
      ELSE
         WRITE(MYUNIT,'(A,I3,A,I10,A)') 'Starting ',NPAR,' parallel MC runs of ',NSTEPS,' steps'
      ENDIF
      WRITE(MYUNIT,'(A,F15.8,A)') 'Temperature will be multiplied by ',SCALEFAC,' at every step'

      IF (RATIOT) THEN
         NULLMOVES=0
         NULLMOVEST=0
      ENDIF

      J1=0
      COMPLETE=.FALSE.
      IF (NSTEPS.EQ.0) COMPLETE=.TRUE. ! otherwise we always take a step


!  Main basin-hopping loop 
      DO WHILE (.NOT.COMPLETE)
         J1=J1+1
         IF (J1.GE.NSTEPS) COMPLETE=.TRUE.
         IF (RELAXRIGIDT) THEN
            IF( (MOD(J1,NRELAXRIGIDR + NRELAXRIGIDA) ==  1) ) THEN
               RIGIDINIT = .TRUE.
               WRITE (MYUNIT, *) 'calling rigid quench'
               CALL CHECKSITES(.TRUE.,COORDS(:,1))
            ELSE IF ( (MOD(J1,NRELAXRIGIDR + NRELAXRIGIDA) == NRELAXRIGIDR + 1) ) THEN
               RIGIDINIT = .FALSE.
               WRITE (MYUNIT, *) 'calling atom quench'
            ENDIF
         ENDIF
         CALL FLUSH(MYUNIT)

!
!  ********************************* Loop over NPAR parallel runs ******************************
!
#ifdef MPI
#else
         DO JP=1,NPAR
#endif

!       If QUCENTRE is specified, move the centre of coordinates
!       to (0,0,0) before taking the next step (improve this so that you
!       can specify where to move the centre like SETCENTRE?
            IF (QUCENTRET) THEN 
!       David mentioned a possible compiler bug causing problems
!       when you just pass a part of the COORDS array so reading the
!       right bit into TEMPCOORDS first and then back out.
               TEMPCOORDS(1:3*NATOMS)=COORDS(1:3*NATOMS,JP)
               CALL CENTRE2(TEMPCOORDS)
               COORDS(1:3*NATOMS,JP)=TEMPCOORDS(1:3*NATOMS)
            ENDIF
            MCTEMP = TEMP(JP)
23          CONTINUE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! START OF STEP TAKING CALLS!                                                                            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   

            !Check which step to take
            CALL WHICH_MOVE(J1)
            !take grouprotation step
            IF (GROUPROTT.AND.DOGROUPROT) CALL GROUPROTSTEP(JP)   
            !rigid body steps
            IF (RIGIDINIT.AND.ROTATERIGIDT.AND.DOROTATERIGID) THEN
               CALL GENRIGID_ROTATE(COORDS(:,JP), ROTATEFACTOR)
            ENDIF
            IF (RIGIDINIT.AND.TRANSLATERIGIDT.AND.DOTRANSLATERIGID) THEN
               CALL GENRIGID_TRANSLATE(COORDS(:,JP),TRANSLATEFACTOR)
            ENDIF
            !hinge moves
            IF (BPHINGET.AND.DOHINGE) THEN 
               ! call potential to get correct base paring (if previous step was rejected)
               CALL UPDATE_POT(COORDS(:,JP))             
               NOBPT = .FALSE. !set control variable to identify problems
               CALL ANALYSE_BP(LOOSEFT, LOOSETT, NUCF, NUCT)
               IF (LOOSEFT.AND.LOOSETT) THEN
                  WRITE(MYUNIT,'(A)') ' mc> Both ends are loose, attempt to move both'
                  CALL MOVE_BOTH_TAILS(COORDS(:,JP), NUCF, NUCT)
               ELSE
                  IF (LOOSEFT) THEN
                     WRITE(MYUNIT,'(A)') ' mc> 5-end is a loose tail, attempt to move it'
                     CALL MOVE_FTAIL(COORDS(:,JP), NATOMS, NUCF, .TRUE.)
                  ELSE IF (LOOSETT) THEN
                     WRITE(MYUNIT,'(A)') ' mc> 3-end is a loose tail, attempt to move it'
                     CALL MOVE_TTAIL(COORDS(:,JP), NATOMS, NUCT, .TRUE.)
                  ELSE
                     !if we don't have base pairs, do a different set of moves
                     IF (NOBPT) THEN
                        WRITE(MYUNIT,'(A)') ' mc> No base pairs found - attempt to move tails towards each other'
                        CALL MOVE_TAILS(COORDS(:,JP),ITERATIONS)
                     ELSE
                        WRITE(MYUNIT,'(A)') ' mc> No loose tail detected, no move attempted, Cartesian step instead'
                        DOCARTSTEP=.TRUE.
                     ENDIF
                  END IF
               END IF
            END IF
            !harmonic constraint move
            IF (HARMONICMOVET.AND.DOHARMONIC) THEN
               ! call potential to get correct base paring (if previous step was rejected)
               CALL UPDATE_POT(COORDS(:,JP))   
               CALL HARMONIC_MOVE(COORDS(:,JP),MOVEDT)
               IF (.NOT.MOVEDT) THEN
                  DOCARTSTEP=.TRUE.
                  WRITE(MYUNIT,*) " mc> No harmonic spring move attempted, do cartesian step instead"
               ENDIF
            ENDIF
            !pulling moves
            IF (PULLMOVET.AND.DOPULL) THEN
               CALL PULL_MOVE(COORDS(:,JP))
               TWISTORPULL=.FALSE.
            ENDIF
            !twisting moves
            IF (TWISTMOVET.AND.DOTWIST) THEN
               IF (TWISTT) THEN
                  WRITE(MYUNIT,*) " mc> twist moves CANNOT be used with a twist potential, deactivating twist moves"
                  TWISTMOVET=.FALSE.
                  DOCARTSTEP=.TRUE.
               ELSE
                  CALL TWIST_MOVE(COORDS(:,JP))
                  TWISTORPULL=.TRUE.
               ENDIF
            ENDIF 
            ! kr366> Gradient modification steps
            IF (STOCHFORCET.AND.DOGRADMODSTEP) THEN
               WRITE(MYUNIT,'(A)') " mc> Use stochastic forces (gradient modification)"
               CALL GRADMOD_STEP(COORDS(:,JP))
            END IF

            ! k2262470> SAXS steps
            IF (SAXSSTEPST.AND.DOSAXSSTEP) THEN
               WRITE(MYUNIT,'(A)') " mc> Attempt SAXS force step"
               CALL CPU_TIME(TSAXS1)
               CALL HIRE_SAXS_FORCE(3*NATOMS,COORDS(:,JP),ESAXS,SAXSFORCE,.TRUE.)
               CALL CPU_TIME(TSAXS2)
               WRITE(MYUNIT,*) " mc> SAXS force: ", DSQRT(SUM(SAXSFORCE(1:3*NATOMS)**2)/(3*NATOMS))
               SAXSSTEPSIZE = RMSLIMITSAXS/MAX(DSQRT(SUM(SAXSFORCE(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100)
               WRITE(MYUNIT,*) " mc> SAXS step size: ", SAXSSTEPSIZE
               WRITE(MYUNIT,*) " mc> SAXS step time: ", TSAXS2-TSAXS1
               COORDS(:,JP) = COORDS(:,JP) + SAXSSTEPSIZE*SAXSFORCE(1:3*NATOMS)
            END IF
            !do Cartesian steps - only if no ther move is attempted!
            IF (DOCARTSTEP) CALL CARTESIAN_SPHERE(COORDS(:,JP), STEP(JP))

            ! Restore atom coordinates if atom is FROZEN or DONTMOVE
            IF(FREEZE) THEN
               DO J2=1,NATOMS
                  IF (FROZEN(J2)) THEN
                     COORDS(3*(J2-1)+1:3*(J2-1)+3,JP)=SAVECOORDS(3*(J2-1)+1:3*(J2-1)+3)
                  ENDIF
               ENDDO
            ENDIF
            ! Redefine rigid bodies if they are being changed by steps
            IF (RIGIDINIT.AND.UPDATERIGIDREFT) THEN
               CALL GENRIGID_UPDATE_REFERENCE(COORDS(:,JP))
            ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! END OF STEP TAKING CALLS!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

            NQ(JP)=NQ(JP)+1
            CALL QUENCH(.FALSE.,JP,ITERATIONS,TIME,QDONE,SCREENC)  
            NQTOT=NQTOT+1
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I0.2,A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') '[',JP, &
     &                       ']Qu ',NQ(JP),' E=',POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
            ELSE
               WRITE(MYUNIT,'(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F11.1)') 'Qu ',NQ(JP),' E=', &
     &                       POTEL,' steps=',ITERATIONS,' RMS=',RMS,' Markov E=',EPREV(JP),' t=',TIME-TSTART
            ENDIF
            CALL FLUSH(MYUNIT)
               
            !lm759> Saving H-bonds in OPEP+HiRE  
            IF(RNAHBT) CALL SAVE_HB(J1,JP)
        
!     TRACKDATA keyword prints the quench energy, markov energy
!     and energy of the current lowest minimum to files for viewing during a run. 

            IF (TRACKDATAT) THEN
               WRITE(MYEUNIT,'(I10,F20.10)') J1,POTEL
               WRITE(MYMUNIT,'(I10,G20.10)') J1,EPREV(JP)
               WRITE(MYBUNIT,'(I10,G20.10)') J1,QMIN(1)
               CALL FLUSH(MYEUNIT)
               CALL FLUSH(MYMUNIT)
               CALL FLUSH(MYBUNIT)
            ENDIF

            ! A series of tests start here to check if a structure should
            ! be allowed into the markov chain. If it fails a test, the ATEST
            ! variable will be set to .FALSE. 
            ATEST=.TRUE.
            IF (COLDFUSION) THEN
               ATEST=.FALSE.
            ENDIF
            COLDFUSION=.FALSE.
            IF (GRADPROBLEMT) THEN
               ATEST=.FALSE.
            ENDIF

            !  Sanity check to make sure the Markov energy agrees with COORDSO. 
            !  Stop if not true.
            IF (DEBUG.AND.(.NOT.FEBHT)) THEN
               CALL POTENTIAL(COORDSO(:,JP),GRAD,OPOTEL,.FALSE.,.FALSE., .FALSE.)
               IF (ABS(OPOTEL-EPREV(JP)).GT.ABS(ECONV)) THEN
                  WRITE(MYUNIT,'(2(A,G20.10))') 'mc> ERROR - energy for coordinates in COORDSO=',OPOTEL, &
     &                                                  ' but Markov energy=',EPREV(JP)
               ELSE
                  WRITE(MYUNIT,'(2(A,G20.10))') 'mc> energy for coordinates in COORDSO=',OPOTEL, &
     &                                                  ' and Markov energy=',EPREV(JP)
               ENDIF
            ENDIF 


            !Accepting criterion
            ! For SAXS forces we use a different comparison
            ! lm59> Only here turn on SAXS force (for OPEP-HIRE interface)
            IF (ATEST.AND.SAXST) THEN
               CALL ADD_SAXS_FORCE(J1,JP,ATEST,MCTEMP)
            ELSE IF (ATEST.AND..NOT.RATIOT) THEN
               CALL CANONICAL_ACC(POTEL,EPREV(JP),ATEST,MCTEMP)
            ENDIF


            IF (RATIOT) THEN
               CALL NULLMOVE(SCREENC,COORDSO(:,JP),ATEST,NULLMOVES,JP)
            ENDIF

            ! Accept or reject step. If the quench did not converge then allow a
            ! potential move, but count it as a rejection in terms of NSUCCESS and
            ! NFAIL. This way we will accept a lower minimum if found, but the steps won;t become so big.
            IF (ATEST) THEN
               IF (DEBUG) THEN
                  WRITE(MYUNIT,334) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
334               FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' ACC')
               ENDIF

               IF (QDONE.EQ.1) THEN
                  NSUCCESS(JP)=NSUCCESS(JP)+1
               ELSE
                  NFAIL(JP)=NFAIL(JP)+1
               ENDIF
               EPREV(JP)=POTEL
               REJSTREAK=0
               COORDSO(1:3*NATOMS,JP)=COORDS(1:3*NATOMS,JP)
            ELSE
               NFAIL(JP)=NFAIL(JP)+1
               REJSTREAK=REJSTREAK+1
               COORDS(1:3*NATOMS,JP)=COORDSO(1:3*NATOMS,JP)
               WRITE(MYUNIT, *) ' mc> Move rejected'
               IF (DEBUG) THEN
                  WRITE(MYUNIT,36) JP,RANDOM,POTEL,EPREV(JP),NSUCCESS(JP),NFAIL(JP)
36                FORMAT('JP,RAN,POTEL,EPREV,NSUC,NFAIL=',I2,3F15.7,2I6,' REJ')
               ENDIF
            ENDIF

            ! Check the acceptance ratio. 
            IF ((MOD(J1,NACCEPT).EQ.0)) THEN
               IF (RATIOT) THEN
                  CALL FIXRATIO(NQ(JP),NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               ELSE
                  CALL ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
               ENDIF
            ENDIF

            TEMP(JP)=TEMP(JP)*SCALEFAC
            
            ! If the move was pulling or twisting, adapt force
            IF ((DOTWIST.AND.TADAPTFT).OR.(DOPULL.AND.PADAPTFT)) THEN
               CALL ADAPT_FORCE(DOTWIST,DOPULL,ATEST)
            ENDIF
               
#ifdef MPI
#else
            IF (EPREV(JP)<EBEST(JP)) EBEST(JP)=EPREV(JP)
#endif
            IF (DUMPINT.GT.0) THEN
               IF (MOD(J1,DUMPINT).EQ.0) THEN
                  CALL DUMPSTATE(NQ(JP),EBEST,BESTCOORDS,JBEST,JP)                  
               ENDIF
            ENDIF

            IF (DUMP_MARKOV) THEN
               ! Dump current Markov state to XYZ for post-processing 
               ! (e.g. characterisation of equilibrium structure)
               J2=J1-1-DUMP_MARKOV_NWAIT
               IF (J2.GE.0 .AND. MOD(J2,DUMP_MARKOV_NFREQ).EQ.0) THEN
                  WRITE(DUMPXYZUNIT(JP),'(I4)') NATOMS
                  WRITE(DUMPXYZUNIT(JP),'(A,I9,A,F20.10)') 'MC step ',J1,' Markov energy ', EPREV(JP)
                  DO J2=1,NATOMS
                     WRITE(DUMPXYZUNIT(JP),'(A,3(1X,F10.5))') 'X ', (COORDSO(3*(J2-1)+J3,JP), J3=1,3)
                  ENDDO
               ENDIF
            ENDIF   
         
#ifdef MPI
            CALL BHPT_EXCHANGE(BETA,J1,NTOT,NEACCEPT)
 
            IF (NQ(1).GT.NSTEPS) EXIT
#else
            IF (NQ(JP).GT.NSTEPS) GOTO 37
         ENDDO
#endif
!  ****************************** End of loop over NPAR parallel runs *****************************
#ifdef MPI
         CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
#endif
         CALL FLUSH(MYUNIT)
      ENDDO                      !end of outer BH loop

37    CONTINUE
#ifdef MPI
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      IF (MYNODE.EQ.0) THEN
         MPI_FINT=.TRUE.
      ELSE
         MPI_FINT=.FALSE.
      ENDIF
      DO
         CALL MPI_BCAST(MPI_FINT,1,MPI_LOGICAL,0,MPI_COMM_WORLD,MPIERR) 
         IF (MPI_FINT) THEN
            WRITE(MYUNIT,'(A)') 'mc> MC steps using MPI finished.'
            GOTO 38
         ENDIF
         WAIT(10)
      ENDDO
38    CONTINUE 
#endif


#ifdef MPI
      CALL BHPT_IO(TIME-TSTART,NQTOT-1,NPCALL-1,NTOT,NEACCEPT)
      WRITE(MYUNIT,10) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),STEP(JP),ASTEP(JP),TEMP(JP)
10    FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)

#else
      DO JP=1,NPAR
         IF (NPAR.GT.1) THEN
            IF (RATIOT) THEN 
               IF (MOD(NQ(JP)-1,NACCEPT).NE.0) CALL FIXRATIO(NQ(JP),NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               CALL FINALRATIO(NULLMOVEST,NSUCCESST,JP,SUM(NQ(:)))
            ENDIF 
            WRITE(MYUNIT,20) JP,NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),STEP(JP),ASTEP(JP),TEMP(JP)
20          FORMAT('[',I2,']Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ELSE
            IF (RATIOT) THEN
               IF (MOD(J1-1,NACCEPT).NE.0) CALL FIXRATIO(J1,NULLMOVES,NSUCCESS,JP,NULLMOVEST,NSUCCESST)
               CALL FINALRATIO(NULLMOVEST,NSUCCESST,JP,SUM(NQ(:)))
            ENDIF
            WRITE(MYUNIT,21) NSUCCESST(JP)*1.0D0/MAX(1.0D0,1.0D0*(NSUCCESST(JP)+NFAILT(JP))),STEP(JP),ASTEP(JP),TEMP(JP)
21          FORMAT('Acceptance ratio for run=',F12.5,' Step=',F12.5,' Angular step factor=',F12.5,' T=',F12.5)
         ENDIF
      ENDDO
#endif
      ! deallocate arrays at the end
      IF (ALLOCATED(TMOVE)) DEALLOCATE(TMOVE)
      IF (ALLOCATED(OMOVE)) DEALLOCATE(OMOVE)
      IF (ALLOCATED(SCREENC)) DEALLOCATE(SCREENC)
      IF (ALLOCATED(SAVECOORDS)) DEALLOCATE(SAVECOORDS)  
      IF (ALLOCATED(TEMPCOORDS)) DEALLOCATE(TEMPCOORDS)
      IF (ALLOCATED(BESTCOORDS)) DEALLOCATE(BESTCOORDS)     
      RETURN
      END SUBROUTINE MC


      SUBROUTINE ACCREJ(NSUCCESS,NFAIL,JP,NSUCCESST,NFAILT)
      USE COMMONS
      USE PREC
      USE grouprotmod
      IMPLICIT NONE
      INTEGER NSUCCESS(NPAR), NFAIL(NPAR), JP, NFAILT(NPAR), NSUCCESST(NPAR)
      REAL(KIND = REAL64) P0,FAC

      ! Calculate the fraction of steps that have succeeded, the acceptance ratio
      P0=1.D0*NSUCCESS(JP)/(1.D0*(NSUCCESS(JP)+NFAIL(JP)))
      ! Set the scaling factor to be used to scale STEP/TEMPERATURE
      ! Calculate the scaling factor using the 'acceptance-ratio method' (ARM)
      IF(ARMT) THEN
         FAC=LOG(ARMA*ACCRAT(JP)+ARMB)/LOG(ARMA*P0+ARMB)
      ELSE
         ! Otherwise just set it to 1.05 
         FAC=1.05D0
      ENDIF
      ! CASE 1: P0 > ACCRAT
      ! If the current acceptance ratio (P0) is larger than the target acceptance ratio (ACCRAT)
      ! we want to increase the step size or reduce the temperature to try to reduce it
      IF (P0.GT.ACCRAT(JP)) THEN
         ! If both the STEP and TEMPERATURE are fixed with FIXBOTH, do nothing
         IF (FIXBOTH(JP)) THEN
         ! Else, if only the STEP is fixed, scale down the TEMPERATURE
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) THEN
               TEMP(JP)=TEMP(JP)/FAC
            ENDIF
         ELSE
            STEP(JP)=STEP(JP)*FAC
            ! GROUPROTATION scaling
            IF(GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) THEN 
               CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,FAC)
            ENDIF
            ! sf344> scale the GENRIGID translation and rotation step sizes
            IF(TRANSLATERIGIDT) TRANSLATEFACTOR=TRANSLATEFACTOR*FAC
            IF(ROTATERIGIDT) ROTATEFACTOR=ROTATEFACTOR*FAC
            ! ASTEP is the maximum angular step size? Not 100% sure, may be a threshold for Morse/LJ...
            ASTEP(JP)=ASTEP(JP)*FAC
         ENDIF
      ! CASE 2: P0 < ACCRAT
      ! If the current acceptance ratio (P0) is smaller than the target acceptance ratio (ACCRAT)
      ! we want to decrease the step size or increase the temperature to try to increase it
      ELSE
         ! If both the STEP and TEMPERATURE are fixed with FIXBOTH, do nothing
         IF (FIXBOTH(JP)) THEN
         ! Else, if only the TEMPERATURE is fixed, scale up the TEMPERATURE
         ELSE IF (FIXSTEP(JP)) THEN
            IF (.NOT.FIXTEMP(JP)) THEN 
               TEMP(JP)=TEMP(JP)*FAC
            ENDIF
         ELSE
            ! This is the usual case for most GMIN runs
            STEP(JP)=STEP(JP)/FAC
            ! GROUPROTATION scaling
            IF (GROUPROTT.AND.(GR_SCALEROT.OR.GR_SCALEPROB)) THEN 
               CALL GROUPROTSCALE(GR_SCALEPROB,GR_SCALEROT,1.0D0/FAC)
            ENDIF
            ! scale the GENRIGID translation and rotation step sizes
            IF(TRANSLATERIGIDT) TRANSLATEFACTOR=TRANSLATEFACTOR/FAC
            IF(ROTATERIGIDT) ROTATEFACTOR=ROTATEFACTOR/FAC           
            ! ASTEP is the maximum angular step size? Not 100% sure, may be a threshold for Morse/LJ...
            ASTEP(JP)=ASTEP(JP)/FAC
         ENDIF
      ENDIF

! Prevent steps from growing out of bounds. The value of 1000 seems sensible, until
! we do something with such huge dimensions?!
      STEP(JP)=MIN(STEP(JP),1.0D3)
      ASTEP(JP)=MIN(ASTEP(JP),1.0D3)
      ! Do some printing to give the user feedback
      IF (MOD(NQ(JP),PRTFRQ).EQ.0) THEN
         IF (NPAR.GT.1) THEN
            WRITE(MYUNIT,'(A,I0.2,A,I6,A,F8.4,A,F8.4)') '[',JP,']Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
         ELSE
            WRITE(MYUNIT,'(A,I6,A,F8.4,A,F8.4)') 'Acceptance ratio for previous ',NACCEPT,' steps=',P0,'  FAC=',FAC
         ENDIF
         IF (FIXBOTH(JP)) THEN
         ELSE IF (FIXSTEP(JP)) THEN
            IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F12.4)') 'Temperature is now:',TEMP(JP)
         ELSE
            IF (NPAR.GT.1) THEN
               WRITE(MYUNIT,'(A,I0.2,A)',ADVANCE='NO') '[',JP,']Steps are now:'
            ELSE
               WRITE(MYUNIT,'(A)',ADVANCE='NO') 'Steps are now:'
            ENDIF
            WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO') '  STEP=',STEP(JP)    
            IF(ASTEP(JP).GT.0.D0) WRITE(MYUNIT,'(A,F10.4)',ADVANCE='NO')'  ASTEP=',ASTEP(JP) 
            IF(TRANSLATERIGIDT) WRITE(MYUNIT,'(A,F10.4)')'  TRANSLATEFACTOR=',TRANSLATEFACTOR
            IF(ROTATERIGIDT) WRITE(MYUNIT,'(A,F10.4)')'  ROTATEFACTOR=',ROTATEFACTOR
            IF(.NOT.FIXTEMP(JP)) WRITE(MYUNIT,'(A,F10.4)') ' Temperature is now:',TEMP(JP)
         ENDIF
      ENDIF

      NSUCCESST(JP)=NSUCCESST(JP)+NSUCCESS(JP)
      NFAILT(JP)=NFAILT(JP)+NFAIL(JP)
      NSUCCESS(JP)=0
      NFAIL(JP)=0 

      RETURN
      END SUBROUTINE ACCREJ

      SUBROUTINE CANONICAL_ACC(ENEW,EOLD,ATEST,MCTEMP)
      IMPLICIT NONE
      REAL(KIND = REAL64), INTENT(IN) :: ENEW, EOLD, MCTEMP
      LOGICAL, INTENT(INOUT) :: ATEST
      REAL(KIND = REAL64) :: DPRAND, RANDOM
    
!  Standard canonical sampling.
      IF (ENEW.LT.EOLD) THEN
         RANDOM=0.0D0
         ATEST=.TRUE.
      ELSE
         RANDOM=DPRAND()
         IF (DEXP(-(ENEW-EOLD)/MAX(MCTEMP,1.0D-100)).GT.RANDOM) THEN
            ATEST=.TRUE.
         ELSE
            ATEST=.FALSE.
         ENDIF
      ENDIF

      RETURN 
      END SUBROUTINE CANONICAL_ACC

      SUBROUTINE WHICH_MOVE(J1)
      USE COMMONS
      USE STOCH_FORCE_STEPS, ONLY: GRADMODFREQ, GRADMODOFFSET, DOGRADMODSTEP, STOCHFORCET
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: J1

      DOGROUPROT = .FALSE.
      DOROTATERIGID = .FALSE.
      DOTRANSLATERIGID = .FALSE.
      DOHINGE = .FALSE.
      DOHARMONIC= .FALSE.
      DOTWIST = .FALSE.
      DOPULL = .FALSE.
      DOCARTSTEP = .FALSE.
     DOGRADMODSTEP = .FALSE.
      !largest move have highest priority
      IF (BPHINGET.AND.MOD(J1,BPHINGEFREQ).EQ.0) THEN
         DOHINGE = .TRUE.
         RETURN
      ENDIF
      !next highest priority for rigid body moves
      IF (ROTATERIGIDT.AND.MOD((J1-ROTRIGIDOFF),ROTATERIGIDFREQ).EQ.0) DOROTATERIGID=.TRUE.
      IF (TRANSLATERIGIDT.AND.MOD((J1-TRANSRIGIDOFF),TRANSLATERIGIDFREQ).EQ.0) DOTRANSLATERIGID=.TRUE.
      IF (DOROTATERIGID.OR.DOTRANSLATERIGID) RETURN
      !next harmonic moves
      IF (HARMONICMOVET.AND.MOD(J1,HARMOVEFREQ).EQ.0) THEN
         DOHARMONIC=.TRUE.
         RETURN
      ENDIF
      !then pulling and twisting moves
      IF (PULLMOVET.AND.MOD((J1-PULLMOFF),PULLMFREQ).EQ.0) DOPULL=.TRUE.
      IF (TWISTMOVET.AND.MOD((J1-TWISTMOFF),TWISTMFREQ).EQ.0) DOTWIST=.TRUE.
      !but as they are big changes, only allow one of them
      IF (DOPULL.AND.DOTWIST) THEN
         IF (TWISTORPULL) THEN !True if last move was twist
            DOTWIST=.FALSE.
         ELSE
            DOPULL=.FALSE.
         ENDIF
         RETURN
      ELSE IF (DOPULL.OR.DOTWIST) THEN
         RETURN
      ENDIF
      !now move to group rotations
      IF (GROUPROTT.AND.MOD(J1,GROUPROTFREQ).EQ.0) THEN
         DOGROUPROT=.TRUE.
         RETURN
      ENDIF
      DOCARTSTEP = .TRUE.

      ! gradient modification steps
      IF (STOCHFORCET.AND.(MOD(J1-GRADMODOFFSET,GRADMODFREQ).EQ.0)) THEN
         DOGRADMODSTEP=.TRUE.
         DOCARTSTEP=.FALSE.
      END IF

      ! SAXS force steps
      IF (SAXSSTEPST.AND.(MOD(J1,SAXSFORCESTEPFREQ).EQ.0)) THEN
         DOSAXSSTEP=.TRUE.
         DOCARTSTEP=.FALSE.
      END IF     
      RETURN
      END SUBROUTINE WHICH_MOVE

      SUBROUTINE SAVE_HB(J1, JP)
        USE HIRE_INTERFACE, ONLY: CALL_HBDAT
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: J1, JP
        REAL(KIND = REAL64) :: HBX(3*NATOMS)

        HBSAVET=(MOD(J1,HBNSTEPS).EQ.0)
        IF(HBSAVET)THEN
            HBX(:)=COORDSO(:,JP)
            CALL CALL_HBDAT(3*NATOMS,HBX)
            HBSAVET=.FALSE.
        ENDIF
      END SUBROUTINE SAVE_HB


     SUBROUTINE ADD_SAXS_FORCE(J1,JP,ATEST,MCTEMP)
        USE COMMONS
        USE HIRE_INTERFACE, ONLY: MODULATE_SAXS
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: J1, JP
        LOGICAL, INTENT(OUT) :: ATEST
        REAL(KIND = REAL64), INTENT(IN) :: MCTEMP
        REAL(KIND = REAL64) :: SAXSX(3*NATOMS), SAXSF(3*NATOMS), SAXSE, SAXSEO
        LOGICAL :: SAXSFORCET
        
        ! potential curve modulation
        SAXSSAVET=(MOD(J1,SAXSNPRINT).EQ.0)
        SAXSFORCET=.TRUE.

        !IF (SAXSMODULT.AND.(J1.LT.SAXSWAVE)) THEN
        IF (SAXSMODULT) THEN
           !SAXSMODI = J1 / SAXSWAVE
           SAXSMODI =  MOD(INT(J1/SAXSNSTEPS),INT(SAXSWAVE))/ SAXSWAVE
           SAXSOFFI = J1 / FLOAT(SAXSNSTEPS)
           !SAXSOFFI = 0.0
           !IF (MOD(J1,SAXSNSTEPS).EQ.0) THEN
           !  SAXSOFFI = 1.0
           !ENDIF
           SAXSFORCET=(MOD(J1,SAXSNSTEPS).EQ.0)
        ELSEIF (.NOT.SAXSMODULT) THEN 
           SAXSFORCET=(MOD(J1,SAXSNSTEPS).EQ.0)
        ELSE
           SAXSFORCET=.FALSE.
           SAXSSAVET=.FALSE.
        ENDIF
        !pass modulation
        CALL MODULATE_SAXS(SAXSINVSIG,SAXSWAVE,SAXSOFFI,SAXSMODI)
        SAXSX(:)=COORDSO(:,JP)
        CALL POTENTIAL(SAXSX, SAXSF, SAXSEO, .FALSE., .FALSE., SAXSFORCET)
        SAXSX(:)=COORDS(:,JP)
        CALL POTENTIAL(SAXSX, SAXSF, SAXSE, .FALSE., .FALSE., SAXSFORCET)
!        SAXSFORCET=.FALSE. 
        SAXSSAVET=.FALSE. 
        CALL CANONICAL_ACC(SAXSE,SAXSEO,ATEST,MCTEMP) !MC here! 
        IF (TRACKDATAT) THEN
           WRITE (MYSEUNIT,'(I10,F20.10)') J1,SAXSE 
!            IF (SAXSMODULT) THEN
!               SAXSFORCET=.TRUE.
!               SAXSMODI = 0.0
!               SAXSOFFI = 1.0
!               CALL MODULATE_SAXS(SAXSINVSIG,SAXSWAVE,SAXSOFFI,SAXSMODI)              
!               CALL POTENTIAL(SAXSX, SAXSF, SAXSE, .FALSE., .FALSE., SAXSFORCET)
!               SAXSFORCET=.FALSE.
!            ENDIF
           IF (ATEST .AND. SAXSFORCET) THEN
            WRITE(MYSMUNIT,'(I10,2F20.10)') J1,SAXSEO,SAXSE 
           ENDIF
           CALL FLUSH(MYSEUNIT)
           CALL FLUSH(MYSMUNIT)
           SAXSFORCET=.FALSE. 
        ENDIF
        RETURN
     END SUBROUTINE ADD_SAXS_FORCE
END MODULE MCmod
