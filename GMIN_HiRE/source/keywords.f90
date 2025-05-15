
      SUBROUTINE KEYWORD

      use COMMONS
      use TWIST_MOD
      use GENRIGID
      USE GROUPROTMOD, ONLY : GROUPROT_INIT
      USE INPUTmod, ONLY : NITEMS, INPUT, READA, READI, READF, REPORT
      USE PORFUNCS

      USE HIRE_INTERFACE, ONLY : SETUP_SAXS, HIRE_HB_INIT
      ! stochastic force steps
      USE STOCH_FORCE_STEPS, ONLY: STOCHFORCET, MULTIPLIER, NMAXMULTSTEPS, NCOMPONENTS, RANDOMFACTORST, &
                                   RMSCUTOFF, JUMPMULTIPLIER, MINMULT, LOWERMOD, RAISEMOD, &
				   GRADMODOFFSET, GRADMODFREQ
      IMPLICIT NONE

      INTEGER IX, I, J1, JP, NPCOUNT, NDUMMY
      INTEGER DATA_UNIT, FUNIT, XUNIT

      LOGICAL YESNO, ENDT

      REAL(KIND = REAL64), ALLOCATABLE :: COORDS1(:)

      CHARACTER(LEN=20) :: UNSTRING
      CHARACTER(LEN=16) :: WORD
      CHARACTER(LEN=10) :: WORD2
      INTEGER :: GROUPOFFSET = 0
      INTEGER :: GETUNIT

      UNSTRING='UNDEFINED'
      NPCOUNT=0

!     >    COORDS and LABELS allocated here
      ALLOCATE(FIXSTEP(1),FIXTEMP(1),FIXBOTH(1),TEMP(1),ACCRAT(1),STEP(1),ASTEP(1),NQ(1),EPREV(1), &
     &     COORDS(3*NATOMSALLOC,1),COORDSO(3*NATOMSALLOC,1))
      DO JP=1,1
         EPREV(JP)=1.0D100
         FIXSTEP(JP)=.FALSE.
         FIXTEMP(JP)=.FALSE.
         FIXBOTH(JP)=.FALSE.
         TEMP(JP)=0.3D0
         ACCRAT(JP)=0.234D0
         STEP(JP)=0.3D0
         ASTEP(JP)=0.3D0
      ENDDO

      ALLOCATE(FROZEN(NATOMSALLOC))
      ! The FROZENRES array is bigger than needed
      ALLOCATE(FROZENRES(NATOMSALLOC))
      DO J1=1,NATOMSALLOC
         FROZEN(J1)=.FALSE.
         FROZENRES(J1)=.FALSE.
      ENDDO

      ALLOCATE(HARMONICFLIST(NATOMSALLOC))
      ALLOCATE(HARMONICR0(3*NATOMSALLOC))
      DO J1=1,NATOMSALLOC
         HARMONICFLIST(J1)=.FALSE.
      ENDDO
! > Generalised rigid body !
      GENRIGIDT = .FALSE.
      ATOMRIGIDCOORDT = .TRUE.
      RIGIDINIT = .FALSE.
      RELAXFQ = .FALSE.
      RELAXRIGIDT = .FALSE.
      NRELAXRIGIDR = 1000000000
      NRELAXRIGIDA = 1000000000

      AACONVERGENCET = .FALSE.
      
!!!!!!!!!! START READING THE DATA FILE HERE !!!!!!!!!!

      CALL FILE_OPEN('data', DATA_UNIT, .FALSE.)

190   CALL INPUT(ENDT, WORD, DATA_UNIT, .TRUE.)
      IF (ENDT .OR. WORD .EQ. 'STOP') THEN
         IF (NPCOUNT.LT.NPAR) THEN
            DO J1=NPCOUNT+1,NPAR
               STEP(J1)=STEP(1)
               ASTEP(J1)=ASTEP(1)
            ENDDO
         ENDIF
         RETURN
      END IF

      IF (WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'.OR.WORD.EQ."!".OR.WORD.EQ."#") THEN
         GOTO 190

!######################!
!HiRE controls!
!######################!

      ELSE IF (WORD.EQ.'HIRE') THEN
         IF(.NOT.ALLOCATED(COORDS1)) ALLOCATE(COORDS1(3*NATOMS))
         IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
         ALLOCATE(COORDS(3*NATOMS,NPAR))
         XUNIT = GETUNIT()
         OPEN(XUNIT, FILE="start")
         READ(XUNIT, *) (COORDS1(I), I=1,3*NATOMS)
         CLOSE(XUNIT)
         WRITE(MYUNIT,'(A)') " keywords> Read coordinates"
         DO J1=1,NPAR
            COORDS(:,J1) = COORDS1(:)
         END DO

      ELSE IF (WORD.EQ.'HIRETOP') THEN
         CALL READA(TOPNAME)

      ELSE IF (WORD.EQ.'HIRESCALE') THEN
         CALL READA(SCALENAME)

      ELSE IF (WORD.EQ.'TITRATION') THEN
         TITRATION = .TRUE.
         CALL READI(TITMETHOD)         ! which titration method (1 - titration steps, 2 - Fernando's code)
         CALL READI(TITNSTEPS)         ! frequency of tritration steps
   
! Needed for running titration code         
!         sed "s/pH/${pH}/" titration_ff.in > cyl.json
 
      ELSE IF (WORD.EQ.'PH') THEN
         CALL READF(PH)
         WRITE(MYUNIT,'(A,F5.2)') ' keyword> Setting PH to ', PH

      ! lm759> HB detection in HiRE - in development
      ELSE IF (WORD.EQ."RNAHB") THEN
          RNAHBT=.TRUE.
          CALL READI(HBNSTEPS)
          CALL HIRE_HB_INIT()

      ! SAXS code for use with HiRE - in development at the moment
      ELSE IF (WORD.EQ.'SAXSPRINT') THEN
          SAXSPRINT=.TRUE.
          CALL READI(SAXSNPRINT)
      ELSE IF (WORD.EQ.'SAXSMODULATION') THEN
          SAXSMODULT=.TRUE.
          CALL READF(SAXSWAVE)
          CALL READF(SAXSINVSIG)
      ELSE IF (WORD.EQ.'SAXSSOLVENT') THEN
          SAXSSOLT = .TRUE.          
          IF (NITEMS.GT.1) THEN
              REFINET=.TRUE.
              CALL READI(NWATLAY)
              CALL READF(WATRAD)
              CALL READF(WATW)
          ENDIF
      ELSE IF (WORD.EQ.'SAXS') THEN
          SAXST =.TRUE.
          CALL READI(SAXSNSTEPS)
          CALL READF(SAXSMAX)
          CALL SETUP_SAXS(SAXST,SAXSPRINT,SAXSMODULT,SAXSINVSIG,SAXSMAX, &
                          SAXSSOLT,REFINET,WATRAD,NWATLAY)

!
!  Keyword for applied static force.
!
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         WRITE(MYUNIT,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE

!
!   Harmonic constraints
!
      ELSE IF (WORD.EQ.'HARMONICF') THEN
         HARMONICF=.TRUE.
         CALL READF(HARMONICSTR)
         DO J1=1,NITEMS-2
            CALL READI(NDUMMY)
            HARMONICFLIST(NDUMMY)=.TRUE.
         ENDDO

!
!  Twist potential
!
      ELSE IF (WORD.EQ.'TWIST') THEN
         TWISTT=.TRUE.
         INQUIRE(FILE='twistgroups',EXIST=YESNO)
         IF (NITEMS.EQ.7) THEN
            NTWISTGROUPS = 1
            ALLOCATE(TWIST_K(1))
            ALLOCATE(TWIST_THETA0(1))
            ALLOCATE(TWIST_ATOMS(4,1))
            CALL READI(TWIST_ATOMS(1,1))
            CALL READI(TWIST_ATOMS(2,1))
            CALL READI(TWIST_ATOMS(3,1))
            CALL READI(TWIST_ATOMS(4,1))
            CALL READF(TWIST_K(1))
            CALL READF(TWIST_THETA0(1))
         ELSE IF (NITEMS.EQ.1 .AND. YESNO) THEN
            OPEN(UNIT=444,FILE='twistgroups',status='unknown')
            WRITE(MYUNIT,*) 'keyword> Reading in twistgroups..'
            READ(444,*) NTWISTGROUPS
            ALLOCATE(TWIST_K(NTWISTGROUPS))
            ALLOCATE(TWIST_THETA0(NTWISTGROUPS))
            ALLOCATE(TWIST_ATOMS(4,NTWISTGROUPS))
            DO J1=1, NTWISTGROUPS
               READ(444,*) TWIST_ATOMS(1,J1),TWIST_ATOMS(2,J1), &
     &                     TWIST_ATOMS(3,J1),TWIST_ATOMS(4,J1), &
     &                     TWIST_K(J1),TWIST_THETA0(J1)
            END DO
         ELSE
            WRITE(MYUNIT,*) 'Provide twistgroups file, or specify details of single dihedral in data'
            STOP
         END IF

! stochastic step taking
! kr366> gradient modificatio steps
      ELSE IF (WORD.EQ.'GRADMODSTEPS') THEN
         STOCHFORCET = .TRUE.
         CALL READI(NCOMPONENTS)
         CALL READI(NMAXMULTSTEPS)
         CALL READF(MULTIPLIER)
         WRITE(MYUNIT,'(A,I6,A,I4,A,F8.2)') " keywords> Use gradient modification steps, modifying ", &
                                            NCOMPONENTS, " gradient components, for a maximum of ", &
                                            NMAXMULTSTEPS, " gradient additions per step and an initial modifier of ", -MULTIPLIER

      ELSE IF (WORD.EQ.'GRADMODADJUST') THEN
         CALL READF(RMSCUTOFF)
         CALL READF(JUMPMULTIPLIER)
         CALL READF(MINMULT)
         WRITE(MYUNIT,'(A,F10.2,A,F8.2,A,F6.2)') " keywords> Use a RMS cutoff of ", RMSCUTOFF, &
                                                 " during gradient modification steps and a multiplier of ", &
                                                 JUMPMULTIPLIER, " to increase multipliers once RMS falls below ", MINMULT
         IF (NITEMS.GT.4) THEN
            CALL READF(LOWERMOD)
            CALL READF(RAISEMOD)
            WRITE(MYUNIT,'(A,F10.2,A,F8.2,A)') " keywords> Multiplier modifications of ", LOWERMOD, " and ",RAISEMOD, &
                                               " to adjust RMS per step"
         END IF

      ELSE IF (WORD.EQ.'GRADMODRANDOM') THEN
         RANDOMFACTORST = .TRUE.

      ELSE IF (WORD.EQ.'GRADMODSTEPFRQ') THEN
	      CALL READI(GRADMODFREQ)
	      CALL READI(GRADMODOFFSET)

! k2262470> SAXS force steps
      ELSE IF (WORD.EQ.'SAXSFORCE') THEN
         SAXSSTEPST = .TRUE.
         CALL READI(SAXSFORCESTEPFREQ)
         CALL READF(RMSLIMITSAXS)
         CALL READF(SAXSMAX)
         CALL SETUP_SAXS(.TRUE.,SAXSPRINT,SAXSMODULT,SAXSINVSIG,SAXSMAX, &
         SAXSSOLT,REFINET,WATRAD,NWATLAY)
!######################!
!BH and minimisation   !
!######################!

      ELSE IF (WORD.EQ.'ACCEPTRATIO') THEN
         IF (NITEMS-1.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of acceptance ratios exceeds NPAR - quit'
            STOP
         ENDIF
         DO J1=1,NITEMS-1
            CALL READF(ACCRAT(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            IF (NPAR.GT.SIZE(ACCRAT)) THEN
               WRITE(MYUNIT,'(A,I10,A,I10)') 'NPAR=',NPAR,' but SIZE(ACCRAT)=',SIZE(ACCRAT)
               WRITE(MYUNIT,'(A,I10,A,I10)') 'Do you need to move the ACCRAT keyword before MPI in data file?'
               STOP
            ENDIF
            DO J1=NITEMS,NPAR
               ACCRAT(J1)=ACCRAT(1)
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'ARM') THEN
         ARMT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ARMA)
         IF (NITEMS.GT.2) CALL READF(ARMB)

!  Interval to check acceptance ratio
      ELSE IF (WORD.EQ.'CHANGEACCEPT') THEN
         CALL READI(NACCEPT)

      ELSE IF (WORD.EQ.'CHECKD') THEN
         CHECKDT = .TRUE.
         IF (NITEMS .GT. 1) CALL READI(CHECKDID)

      ELSE IF (WORD.EQ.'COLDFUSION') THEN
         CALL READF(COLDFUSIONLIMIT)

      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.

!  Initial guess for diagonal matrix elements in LBFGS.
      ELSE IF (WORD.EQ.'DGUESS') THEN
         CALL READF(DGUESS)

      ELSE IF (WORD.EQ.'EDIFF') THEN
         CALL READF(ECONV)

!  Fix temperature and step size
      ELSE IF (WORD.EQ.'FIXBOTH') THEN
         IF (NITEMS.EQ.1) THEN
            FIXBOTH(1)=.TRUE.
            IF (NPAR.GT.1) THEN
               DO J1=2,NPAR
                  FIXBOTH(J1)=.TRUE.
               ENDDO
            ENDIF
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX)
               FIXBOTH(IX)=.TRUE.
            ENDDO
         ENDIF

!  Fix step size
      ELSE IF (WORD.EQ.'FIXSTEP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXSTEP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX)
               FIXSTEP(IX)=.TRUE.
            ENDDO
         ENDIF

!  Fix temperature
      ELSE IF (WORD.EQ.'FIXTEMP') THEN
         IF (NITEMS.EQ.1) THEN
            FIXTEMP(1)=.TRUE.
         ELSE
            DO J1=1,NITEMS-1
               CALL READI(IX)
               FIXTEMP(IX)=.TRUE.
            ENDDO
         ENDIF

      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)

!  Inertia difference criterion - no longer used for distinguishing stationary points!
      ELSE IF (WORD.EQ.'ITOL') THEN
         CALL READF(IDIFFTOL)

      ELSE IF (WORD.EQ.'MAXBFGS') THEN
         CALL READF(MAXBFGS)

      ELSE IF (WORD.EQ.'MAXERISE') THEN
         CALL READF(MAXERISE)
         MAXERISE_SET=.TRUE.

      ELSE IF (WORD.EQ.'MAXIT') THEN
         CALL READI(MAXIT)
         IF (NITEMS.GT.2) THEN
            CALL READI(MAXIT2)
         ENDIF

      ELSE if (WORD.EQ.'QTESTCONV') THEN
         CALL READF(QTESTMAX)

      ELSE IF (WORD.EQ.'RADIUS') THEN
         CALL READF(RADIUS)
         
!  integer seed for random number generator.
      ELSE IF (WORD.EQ.'RANSEED') THEN
         RANSEEDT=.TRUE.
         NDUMMY=0                                   !need to set a default
         CALL READI(NDUMMY)
         CALL SDPRND(NDUMMY+MYNODE)
         CALL SDPRND_UNIVERSAL(NDUMMY+NPAR)
         WRITE(MYUNIT,'(A,I8)') 'keywords> Random seed = ',NDUMMY+MYNODE

      ELSE IF (WORD.EQ.'RATIO') THEN
         RATIOT=.TRUE.
         CALL READF(SRATIO)
         CALL READF(TRATIO)
         FIXBOTH=.TRUE.

!  Restore the state of a previous GMIN run from dumpfile.
      ELSE IF (WORD.EQ.'RESTORE') THEN
         RESTORET=.TRUE.
         CALL READA(DUMPFILE)

      ELSE IF (WORD.EQ.'SLOPPYCONV') THEN
         CALL READF(BQMAX)

!  Sparse Hessian.
      ELSE IF (WORD.EQ.'SPARSE') THEN
         SPARSET=.TRUE.
         IF (NITEMS.GT.1) CALL READF(ZERO_THRESH)


      ELSE IF (WORD.EQ.'STEPS') THEN
         CALL READI(MCSTEPS)
         IF (NITEMS.GT.2) CALL READF(TFAC)
         
!Simulation temperature
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         DO J1=1,NITEMS-1
            CALL READF(TEMP(J1))
         ENDDO
         IF (NITEMS-1.LT.NPAR) THEN
            DO J1=NITEMS,NPAR
               TEMP(J1)=TEMP(1)
            ENDDO
         ENDIF


      ELSE IF (WORD.EQ.'TIGHTCONV') THEN
         CALL READF(CQMAX)

!  Number of BFGS updates before resetting, default=4
      ELSE IF (WORD.EQ.'UPDATES') THEN
         CALL READI(MUPDATE)

!######################!
!Step taking routines  !
!######################!

! Read in the maximum initial step size, factor for determining angular
! moves, and for rigid bodies the angular step size and the size of the
! blocks for Cartesian and angular moves.
!
! For parallel runs different values can be used for different runs by
! adding additional "STEP" lines to the data file. Otherwise the
! parameters for subsequent parallel runs are set to the values for the
! first one.
      ELSE IF (WORD.EQ.'STEP') THEN
         NPCOUNT=NPCOUNT+1
         IF (NPCOUNT.GT.NPAR) THEN
            WRITE(MYUNIT,'(A)') 'Number of STEP lines exceeds NPAR - quit'
            STOP
         ENDIF
         CALL READF(STEP(NPCOUNT))
         CALL READF(ASTEP(NPCOUNT))


! Group rotation moves (now for both AMBER and CHARMM!
      ELSE IF (WORD.EQ.'GROUPROTATION') THEN
         YESNO=.FALSE.
         INQUIRE(FILE='atomgroups',EXIST=YESNO)
         IF (YESNO) THEN
            GROUPROTT=.TRUE.
            WRITE(MYUNIT,'(A)') ' keyword> group rotation moves enabled'
          
         ELSE
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: atom groups must be defined in atomgroups file'
            STOP
         ENDIF
         IF (NITEMS.GT.1) CALL READI(GROUPROTFREQ)
! if the frequency is 0, we need to disable the moves to present a divide by 0!
         IF(GROUPROTFREQ.EQ.0) THEN
            GROUPROTT=.FALSE.
            WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of GROUPROTATION moves set to 0 - moves DISABLED!'
         ENDIF
! Specify GROUPROTATION move scaling mode amd inform the user
         IF (NITEMS.GT.2) THEN
            CALL READA(GR_SCALEMODE)
            IF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEROT') THEN
               GR_SCALEROT=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION amplitudes will be scaled'
            ELSEIF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEPROB') THEN
               GR_SCALEPROB=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION selection probabilities will be scaled'
            ELSEIF(TRIM(ADJUSTL(GR_SCALEMODE)).EQ.'SCALEBOTH') THEN
               GR_SCALEROT=.TRUE.
               GR_SCALEPROB=.TRUE.
               WRITE(MYUNIT,'(A)') ' keyword> GROUPROTATION amplitudes and selection probabilities will be scaled'
            ENDIF
         ENDIF
! Read atom offset for group definitions in atomgroups
         IF (NITEMS.GT.3) CALL READI(GROUPOFFSET)
! Now call the setup subroutine
         CALL GROUPROT_INIT(GROUPOFFSET)  



! Keyword to suppress output of group rotation moves
      ELSE IF (WORD .EQ. 'QUIETGROUPROT') THEN
         GROUPROT_SUPPRESS = .TRUE.

! Keyword to prevent rotation of paired bases
      ELSE IF (WORD .EQ. 'SKIPBPROT') THEN
         SKIPBPT = .TRUE.

! Global hinge moves for loose 3' or 5' tails
      ELSE IF (WORD .EQ. 'BPHINGE') THEN
         BPHINGET = .TRUE.
         CALL READI(BPHINGEFREQ)
         CALL READI(NLOOSE)
         IF (NLOOSE.LT.3) THEN
            WRITE(MYUNIT,*) " keywords> the number of loose nucleotides for a free tail must be at least 3"
            STOP
         END IF
         CALL READF(BPTHRESH)

      ELSE IF (WORD .EQ. 'NOBPSPRING') THEN
         CALL READF(BPDIST)
         CALL READF(BPSTRENGTH)
         CALL READI(NBPHARMOVE)

! Move adding spring to form pull free bases together
      ELSE IF (WORD .EQ. 'HARMONICMOVES') THEN
         HARMONICMOVET = .TRUE.
         CALL READF(HMDIST)
         CALL READF(HMKF)
         CALL READI(HARMOVEFREQ)

! Move adding pulling force
      ELSE IF (WORD .EQ. 'PULLMOVES') THEN
         PULLMOVET = .TRUE.
         CALL READF(PULLMF)
         CALL READI(PULLMFREQ)
         IF (NITEMS.GT.2) THEN
            CALL READI(PULLMOFF)
         ENDIF

! Move adding twisting force
      ELSE IF (WORD .EQ. 'TWISTMOVES') THEN
         TWISTMOVET = .TRUE.
         CALL READF(TWISTMF)
         CALL READI(TWISTMFREQ)
         IF (NITEMS.GT.2) THEN
            CALL READI(TWISTMOFF)
         ENDIF

! Adapt pulling force
      ELSE IF (WORD .EQ. 'ADAPTTWISTF') THEN
         TADAPTFT = .TRUE.                            
         CALL READF(TLOWERF)
         CALL READF(TUPPERF)   
         CALL READF(TADAPTSCALE)

! Adapt pulling force
      ELSE IF (WORD .EQ. 'ADAPTPULLF') THEN
         PADAPTFT = .TRUE.                            
         CALL READF(PLOWERF)
         CALL READF(PUPPERF)   
         CALL READF(PADAPTSCALE)

!######################!
!Rigid bodies+freezing !
!######################!

      ELSE IF (WORD.EQ.'AACONVERGENCE') THEN
         AACONVERGENCET = .TRUE.

      ELSE IF (WORD.EQ.'HYBRIDMIN') THEN
         IF (.NOT.RIGIDINIT) THEN
            WRITE(MYUNIT,'(A)') ' keyword> ERROR: HYBRIDMIN can only be used with RIGIDINIT!'
            STOP
         ENDIF
         CALL READF(EPSRIGID)
         HYBRIDMINT=.TRUE.
         WRITE(MYUNIT,'(A)') ' keyword> Using hybrid rigid body/all-atom minimisation'
         WRITE(MYUNIT,'(A,F20.10)') ' HYBRIDMIN> Rigid body RMS force target= ',EPSRIGID
         WRITE(MYUNIT,'(A)') ' HYBRIDMIN> Final quenches will be done atomistically'

! redefine rigid body sites every NRelaxRigid step   !
      ELSE IF (WORD.EQ.'NRELAXRIGID') THEN
         RELAXRIGIDT = .TRUE.
         CALL READI(NRELAXRIGIDR)
         CALL READI(NRELAXRIGIDA)

! use atom coords during final quench   !
      ELSE IF (WORD.EQ.'RELAXFINALQUENCH') THEN
         RELAXFQ = .TRUE.

! Generalised rigid body   
      ELSE IF (WORD.EQ.'RIGIDINIT') THEN
         RIGIDINIT = .TRUE.
         ATOMRIGIDCOORDT = .TRUE.

! Rigid body rotation moves. Each rigid body is randomly rotated about its COM every ROTATERIGIDFREQ steps.
!        ROTATEFACTOR scales the maximum rotation with 1.0 being complete freedom to rotate.
      ELSE IF (WORD.EQ.'ROTATERIGID') THEN
         ROTATERIGIDT=.TRUE.
! Read ROTATERIGIDFREQ
         IF (NITEMS.GT.1) CALL READI(ROTATERIGIDFREQ)
! Read in ROTATEFACTOR
         IF (NITEMS.GT.2) CALL READF(ROTATEFACTOR)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body rotation moves enabled'
         WRITE(MYUNIT,'(A,I2,A)') ' keyword> Rigid bodies will be rotated every ',ROTATERIGIDFREQ,' steps'
         WRITE(MYUNIT,'(A,F20.10)') ' ROTATERIGID> Rigid body ROTATEFACTOR =',ROTATEFACTOR

! Rigid body translation moves
!        TRANSLATEFACTOR sets the maximum translation distance
      ELSE IF (WORD.EQ.'TRANSLATERIGID') THEN
         TRANSLATERIGIDT=.TRUE.
! Read TRANSLATERIGIDFREQ
         IF (NITEMS.GT.1) CALL READI(TRANSLATERIGIDFREQ)
! Read in TRANSLATEFACTOR
         IF (NITEMS.GT.2) CALL READF(TRANSLATEFACTOR)
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body translation moves enabled'
         WRITE(MYUNIT,'(A,I2,A)') ' keyword> Rigid bodies will be translated every ',TRANSLATERIGIDFREQ,' steps'
         WRITE(MYUNIT,'(A,F20.10)') ' TRANSLATERIGID> Rigid body TRANSLATEFACTOR =',TRANSLATEFACTOR

      ! offset for step taking frequency
      ELSE IF (WORD.EQ.'TRANSRIGIDOFF') THEN
         CALL READI(TRANSRIGIDOFF)

      ! offset for step taking frequency
      ELSE IF (WORD.EQ.'ROTRIGIDOFF') THEN
         CALL READI(ROTRIGIDOFF)


! Update the reference coordinates for the generalised rigid bodies after a step has been taken.
!        This allows steps to be taken WITHIN the rigid bodies, although HYBRIDMIN should also be used
!        as any bad conformation intrduced by these step will otherwise be frozen into the rigid bodies.
      ELSE IF (WORD.EQ.'UPDATERIGIDREF') THEN
         WRITE(MYUNIT,'(A)') ' keyword> Rigid body reference coordinates will be updated after each step'
         WRITE(MYUNIT,'(A)') ' UPDATERIGIDREF> WARNING: make sure HYBRIDMIN is enabled!'
         UPDATERIGIDREFT=.TRUE.


!  Frozen atoms.
      ELSE IF (WORD.EQ.'FREEZE') THEN
         FREEZE=.TRUE.
         IF(NITEMS.GT.1) THEN
            DO J1=1,NITEMS-1
               NFREEZE=NFREEZE+1
               CALL READI(NDUMMY)
               FROZEN(NDUMMY)=.TRUE.
            ENDDO
         ELSE
      ! 2nd input mode: frozen atoms are listed in a file called 'frozen'. The first line
      ! of the file contains the number of frozen atoms, subsequent lines contain atoms to freeze
            INQUIRE(FILE='frozen',EXIST=YESNO)
            IF (YESNO) THEN
                FUNIT=GETUNIT()
                OPEN(FUNIT,FILE='frozen',STATUS='OLD')
                READ(FUNIT,*) NFREEZE
                DO J1=1,NFREEZE
                   READ(FUNIT,*) NDUMMY
                   FROZEN(NDUMMY)=.TRUE.
                ENDDO
            ELSE
                WRITE (*,'(A)') ' ERROR: FREEZE specified incorrectly'
                WRITE(*,*) "Specify frozen atoms either on the keyword line, or in a file called 'frozen'"
                STOP
           ENDIF
         ENDIF

! unfreeze everything at the final quenches
      ELSE IF (WORD.EQ.'UNFREEZEFINALQ') THEN
        UNFREEZEFINALQ=.TRUE.

! Frozen residues (to be converted to frozen atoms)
      ELSE IF (WORD.EQ.'FREEZERES') THEN
         FREEZE=.TRUE.
         FREEZERES=.TRUE.
! The FROZENRES array is then filled with the residue number from the
! data file
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.TRUE.
         ENDDO
! Finally, the frozen residue numbers are converted into frozen atom
! numbers.



! Freezing EVERYTHING and then permitting small parts to move
! This is useful for large system to prevent the data file getting silly
      ELSEIF (WORD.EQ.'FREEZEALL') THEN
         FREEZE=.TRUE.
         FREEZEALL=.TRUE.
         NFREEZE=NATOMSALLOC
         DO J1=1,NATOMSALLOC
            FROZEN(J1)=.TRUE.
            FROZENRES(J1)=.TRUE.
         ENDDO

! Things are then UNFROZEN using the UNFREEZE and UNFREEZERES keywords
! This is only a valid keyword if FREEZEALL is also specified
      ELSEIF ((WORD.EQ.'UNFREEZE').AND.FREEZEALL) THEN
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.FALSE.
            NFREEZE=NFREEZE-1
         ENDDO

      ELSEIF ((WORD.EQ.'UNFREEZERES').AND.FREEZEALL) THEN
         UNFREEZERES=.TRUE.
! Set the right parts of the FROZENRES array to FALSE
         DO J1=1,NITEMS-1
            CALL READI(NDUMMY)
            FROZENRES(NDUMMY)=.FALSE.
         ENDDO
         FREEZERES=.TRUE.


!######################!
!Free energy BH        !
!######################!

     ELSE IF (WORD.EQ.'FEBH') THEN
         CALL READF(FETEMP)
         FEBHT = .TRUE.
         FE_FILE_UNIT = GETUNIT()
         OPEN(UNIT = FE_FILE_UNIT, FILE = 'free_energy', STATUS = 'REPLACE')
         WRITE(FE_FILE_UNIT, '(6A20)') '       Quench       ', '  Potential energy  ', &
     &   '   rot/vib terms    ', '    Free energy     ', '   Markov energy    ', '        Time        '
         WRITE(FE_FILE_UNIT, '(6A20)') ' ------------------ ', ' ------------------ ', &
     &   ' ------------------ ', ' ------------------ ', ' ------------------ ', ' ------------------ '
         IF (NITEMS .GT. 2) THEN
             CALL READA(WORD2)
             WORD2 = TRIM(ADJUSTL(WORD2))
             IF (WORD2 .EQ. 'SPARSE') SPARSET = .TRUE.
         END IF
         IF (NITEMS .GT. 3) THEN
             CALL READF(ZERO_THRESH)
         END IF

! for use with free energy basin-hopping
! this keyword ensures that minimisation is repeated to a stricter convergence criterion if
! the zero and non-zero eigenvalues are not separated by MIN_ZERO_SEP
      ELSE IF (WORD.EQ.'MIN_ZERO_SEP') THEN
         IF (NITEMS.GT.1) THEN
             CALL READF(MIN_ZERO_SEP)
! Refers to frequencies, but checks eigenvalues.
             MIN_ZERO_SEP = MIN_ZERO_SEP**2
         END IF
         IF (NITEMS.GT.2) THEN
             CALL READI(MAX_ATTEMPTS)
         END IF

! Planck constant in prevailing units. Seems necessary for FEBH based calculations.
      ELSE IF (WORD.EQ.'PLANCK') THEN
         CALL READF(PLANCK)

! If USEFRQS is specified we use the quantum partition function.
      ELSE IF(WORD.EQ.'USEFRQS') THEN
         USEFRQS=.TRUE.

! Whether to include rotational partition function for FEBH
      ELSE IF (WORD.EQ.'USEROT') THEN
         USEROT=.TRUE.

!######################!
!Parallel setup        !
!######################!

! PT basin-hopping. This keyword is simply used to read in PTTMIN, PTTMAX, and EXCHPROB.
! It is used in conjunction with MPI to decide if this is a BHPT run in mc.F.
      ELSE IF (WORD.EQ.'BHPT') THEN
         CALL READF(PTTMIN)
         CALL READF(PTTMAX)
         CALL READF(EXCHPROB)
         IF (EXCHPROB.GE.1.D0) THEN
            EXCHINT=INT(EXCHPROB)
            EXCHPROB=1.D0/EXCHPROB
         ELSEIF (EXCHPROB.GT.0.D0) THEN
            EXCHINT=INT(1.D0/EXCHPROB)
         ELSE
            EXCHINT=10000000
         ENDIF
         IF (NITEMS.GT.4) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RANDOM') PTRANDOM=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'INTERVAL') PTINTERVAL=.TRUE.
         ELSE
            PTRANDOM=.TRUE.
         ENDIF
         IF (NITEMS.GT.5) THEN
            CALL READA(UNSTRING)
            WRITE(*,*)UNSTRING
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SINGLE') PTSINGLE=.TRUE.
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SETS') PTSETS=.TRUE.
         ELSE
            PTSINGLE=.TRUE.
         ENDIF
         IF (NITEMS.GT.6) THEN
            CALL READI(NDUMMY)
            CALL SDPRND(NDUMMY)
         ENDIF

!  MPI keyword
      ELSE IF (WORD.EQ.'MPI') THEN
         MPIT=.TRUE.
         DEALLOCATE(FIXSTEP,FIXTEMP,FIXBOTH,TEMP,ACCRAT,STEP,ASTEP,COORDS,NQ,EPREV,COORDSO)

         ALLOCATE(FIXSTEP(NPAR),FIXTEMP(NPAR),FIXBOTH(NPAR),TEMP(NPAR),ACCRAT(NPAR),STEP(NPAR),ASTEP(NPAR), &
     &        COORDS(3*NATOMSALLOC,NPAR),NQ(NPAR),EPREV(NPAR),COORDSO(3*NATOMSALLOC,NPAR))
         DO JP=1,NPAR
            EPREV(JP)=1.0D100
            FIXSTEP(JP)=.FALSE.
            FIXTEMP(JP)=.FALSE.
            FIXBOTH(JP)=.FALSE.
            TEMP(JP)=0.3D0
            ACCRAT(JP)=0.5D0
            STEP(JP)=0.3D0
            ASTEP(JP)=0.3D0
         ENDDO



!######################!
!Coords manipulation   !
!######################!

      ELSE IF (WORD.EQ.'CENTRE') THEN
         CENT=.TRUE.

! SETCENTRE moves the centre of coordinates to the specified
! location before the initial quench is done.
      ELSE IF (WORD.EQ.'SETCENTRE') THEN
         SETCENT=.TRUE.
         IF (NITEMS.EQ.2) THEN
            CALL READF(CENTX)
         ELSE IF (NITEMS.EQ.3) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
         ELSE IF (NITEMS.EQ.4) THEN
            CALL READF(CENTX)
            CALL READF(CENTY)
            CALL READF(CENTZ)
         ENDIF
      
! QUCENTRE moves the centre of coordinates to (0,0,0)
! before each step is taken.
      ELSE IF (WORD.EQ.'QUCENTRE') THEN
         QUCENTRET=.TRUE.


!######################!
!Output controls       !
!######################!

      ELSE IF (WORD.EQ.'DUMPINT') THEN
         CALL READI(DUMPINT)

      ELSE IF (WORD.EQ.'DUMP_MARKOV') THEN
         DUMP_MARKOV=.TRUE.
         CALL READI(DUMP_MARKOV_NWAIT)
         CALL READI(DUMP_MARKOV_NFREQ)

      ELSE IF (WORD.eq.'DUMPMIN') THEN
        DUMPMINT=.TRUE.
        WRITE(MYUNIT,'(A)') ' keywords> The SAVE lowest minima will be dumped every DUMPINT steps as dumpmin.x files'      

      ELSE IF (WORD.EQ.'DUMPPIGS') THEN
        DUMPPIGST = .TRUE.
        IF (NITEMS.GT.1) THEN
           CALL READF(EPIGSLIM)
        ENDIF
        WRITE(MYUNIT,'(A)') ' keywords> Minimisation output for PiGS'

      ELSE IF (WORD.eq.'DUMPSTRUCTURES') THEN
        DUMPSTRUCTURES=.TRUE.
        WRITE(MYUNIT,'(A)') ' keywords> Final structures will be dumped in different formats (.rst, .xyz, .pdb)'

      ELSE IF (WORD.EQ.'PRTFRQ') THEN
         CALL READI(PRTFRQ)

      ELSE IF (WORD.EQ.'SAVE') THEN
         CALL READI(NSAVE)

      ELSE IF (WORD.EQ.'TRACKDATA') THEN
         TRACKDATAT=.TRUE.

      ELSE
         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF
      CALL FLUSH(MYUNIT)

      GOTO 190

      RETURN
      END
