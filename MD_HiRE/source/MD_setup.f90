MODULE MD_SETUP

   CONTAINS
      SUBROUTINE START_TRACKING()
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE MD_COMMONS
         IMPLICIT NONE
         CALL FILE_OPEN("md_energy.log",EUNIT,.TRUE.)
         CALL FILE_OPEN("md_coords.xyz",XUNIT,.TRUE.)
         IF (RMSDT) CALL FILE_OPEN("md_rmsd.log",RUNIT,.TRUE.)
      END SUBROUTINE START_TRACKING

      SUBROUTINE SETUP_POTENTIAL()
         USE MD_COMMONS, ONLY: MYUNIT, NATOMS, NOPT, TOPNAME, SCALEDATNAME, COORDSFILE, &
                               MASSES, COORDS, MININITIAL, RESTARTSIMT, ATNAMES
         USE HIRE_INTERFACE, ONLY: HIRE_INITIALISE, PASS_HIRE_MASSES, PASS_PARTICLE_NAMES
         USE MD_UTILS, ONLY: ALLOC_COMMONS, RUNMIN
         USE FILE_UTILS, ONLY: FILE_EXIST, FILE_OPEN
         IMPLICIT NONE
         INTEGER :: J, XUNIT
         !first initialise the HiRE interface
         CALL HIRE_INITIALISE(TOPNAME, SCALEDATNAME, NATOMS)

         NOPT = 3*NATOMS

         ! allocate the relevant arrays
         CALL ALLOC_COMMONS()


         IF (.NOT.RESTARTSIMT) THEN
            ! get coordinates
            IF (FILE_EXIST(COORDSFILE)) THEN
               CALL FILE_OPEN(COORDSFILE,XUNIT,.FALSE.)
               READ(XUNIT, *) (COORDS(J), J=1,3*NATOMS)
               CLOSE(XUNIT)
               ! minimise coordinates
               IF (MININITIAL) THEN
                  CALL RUNMIN(COORDS)
               END IF
            ELSE
               WRITE(MYUNIT,*) " setup> Cannot locate input file for coordinates - ", COORDSFILE
               STOP
            END IF
         END IF

         ! get particle masses
         CALL PASS_HIRE_MASSES(NATOMS,MASSES)
         !get particle names
         CALL PASS_PARTICLE_NAMES(NATOMS,ATNAMES)
         !get elements
         CALL PARTS2ELS()
      END SUBROUTINE SETUP_POTENTIAL

      SUBROUTINE PARTS2ELS()
         USE MD_COMMONS, ONLY: NATOMS, ATNAMES, ELEMENTS
         IMPLICIT NONE
         INTEGER :: I
         CHARACTER(LEN=4) :: NAME

         DO I=1,NATOMS
            NAME = ATNAMES(I)
            IF ((NAME.EQ."C").OR.(NAME.EQ."P").OR.(NAME.EQ."O")) THEN
               ELEMENTS(I) = TRIM(ADJUSTL(NAME))
            ELSE
               ELEMENTS(I) = "X"
            END IF
         END DO
      END SUBROUTINE PARTS2ELS

      SUBROUTINE READ_SETTINGS()
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE INPUTMOD, ONLY: INPUTKW
         IMPLICIT NONE
         INTEGER :: PARAMUNIT
         LOGICAL :: EOFT
         CHARACTER(25) :: KEYWORD
            
         ! open simdata file    
         CALL FILE_OPEN("mddata", PARAMUNIT, .FALSE.)
            
         !loop over lines in file
         EOFT = .FALSE.
         DO WHILE (.NOT. EOFT)
            CALL INPUTKW(EOFT, KEYWORD, PARAMUNIT, .TRUE.)         
            IF (EOFT) THEN
               EXIT
            ELSE
               CALL SETKEYS(KEYWORD)
            ENDIF
         END DO
         CLOSE(PARAMUNIT)
      END SUBROUTINE READ_SETTINGS

      SUBROUTINE SETKEYS(WORD)
         USE INPUTMOD
         USE MD_COMMONS                  ! global variables
         USE MOD_THERMALISE, ONLY: NTHERMALISE, NEQUIL, NCENTRE, NRMANG, NRESCALE, NEQDUMPE
         USE MINIMISATION, ONLY: ITMAX, MUPDATE, EPS, DUMPINTMIN, DUMPMINCOORDST
         ! USE FILE_UTILS, ONLY: FILE_EXIST
        
         IMPLICIT NONE
         CHARACTER(25), INTENT(IN) :: WORD
         CHARACTER(1) :: CONTINUEDUMMY = "F"
        
         ! Keyword IF clause - first is comments, last is unrecognised command,
         ! everything else in alphabetical order
         ! Labels for each letter is: ! LETTER #A
         ! For each keyword add comment of format:
         ! ! Keyword: COMMENT
         ! ! Added: 11/03/2021 (kr366), last modified: 11/03/2021 (kr366)
         ! ! Description: Keyword to comment lines
         ! Undocumented keywords will be deleted
        
         ! Keyword: COMMENT
         ! Added: 11/03/2021 (kr366), last modified: 11/03/2021 (kr366)
         ! Description: Keyword to comment lines     
         IF (WORD.EQ.'COMMENT'.OR.WORD.EQ."!".OR.WORD.EQ."#") THEN
            RETURN

         !+++++++++++++++!   
         ! LETTER A      !
         !+++++++++++++++!           

         !+++++++++++++++!   
         ! LETTER B      !
         !+++++++++++++++!

         !+++++++++++++++!   
         ! LETTER C      !
         !+++++++++++++++!     
         
         ! Keyword: COORDINATES
         ! Added: 02/11/2022 (kr366), last modified: 02/11/2022 (kr366)
         ! Description: Interval to dump coordinates
         ELSE IF (WORD .EQ. 'COORDINATES') THEN
            CALL READA(COORDSFILE) 

         !+++++++++++++++!   
         ! LETTER D      !
         !+++++++++++++++! 
          
         ! Keyword: DUMPCOORDS
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Interval to dump coordinates
         ELSE IF (WORD .EQ. 'DUMPCOORDS') THEN
            CALL READI(NDUMPX) 

         ! Keyword: DUMPENERGY
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Interval to dump energy
         ELSE IF (WORD .EQ. 'DUMPENERGY') THEN
            CALL READI(NDUMPE) 

         ! Keyword: DUMPPDB
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Interval to dump coordinates as pdb files
         ELSE IF (WORD .EQ. 'DUMPPDB') THEN
            CALL READI(NDUMPP) 
            DUMPPDBT = .TRUE.

         ! Keyword: DUMPRST
         ! Added: 30/11/2022 (kr366), last modified: 30/11/2022 (kr366)
         ! Description: Interval to write restart file
         ELSE IF (WORD .EQ. 'DUMPRST') THEN
            CALL READI(NDUMPRST)            

         !+++++++++++++++!   
         ! LETTER E      !
         !+++++++++++++++!
            
         !+++++++++++++++!   
         ! LETTER F      !
         !+++++++++++++++!
            
         !+++++++++++++++!   
         ! LETTER G      !
         !+++++++++++++++!

         ! Keyword: GAMMA
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Value for friction parameter gamma
         ELSE IF (WORD .EQ. 'GAMMA') THEN
            CALL READF(GAMMA) 

         !+++++++++++++++!   
         ! LETTER H      !
         !+++++++++++++++!            

         !+++++++++++++++!   
         ! LETTER I      !
         !+++++++++++++++!
         
         !+++++++++++++++!   
         ! LETTER J      !
         !+++++++++++++++!     
         
         !+++++++++++++++!   
         ! LETTER K      !
         !+++++++++++++++!     
         
         !+++++++++++++++!   
         ! LETTER L      !
         !+++++++++++++++!      
               
         !+++++++++++++++!   
         ! LETTER M      !
         !+++++++++++++++!

         ! Keyword: MDMODE
         ! Added: 28/11/2022 (kr366), last modified: 28/11/2022 (kr366)
         ! Description: MD method to be used: VV for Velocity Verlet or LD for Langevin Dynamics
         ELSE IF (WORD .EQ. 'MDMODE') THEN
            CALL READA(MDMETHOD)

         ! Keyword: MDSTEPS
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Number of MD steps to be taken
         ELSE IF (WORD .EQ. 'MDSTEPS') THEN
            CALL READI(MDSTEPS) 

         ! Keyword: MININIT
         ! Added: 10/11/2022 (kr366), last modified: 10/11/2022 (kr366)
         ! Description: Run minimisation for initial structure
         ELSE IF (WORD .EQ. 'MININIT') THEN
            MININITIAL = .TRUE.
            CALL READI(ITMAX)
            CALL READI(MUPDATE)
            CALL READF(EPS)
         
         ELSE IF (WORD .EQ. 'MINTRACK') THEN
            DUMPMINCOORDST = .TRUE.
            CALL READI(DUMPINTMIN)

         !+++++++++++++++!   
         ! LETTER N      !
         !+++++++++++++++!

         !+++++++++++++++!   
         ! LETTER O      !
         !+++++++++++++++!    

         !+++++++++++++++!   
         ! LETTER P      !
         !+++++++++++++++!

         !+++++++++++++++!   
         ! LETTER Q      !
         !+++++++++++++++!

         !+++++++++++++++!   
         ! LETTER R      !
         !+++++++++++++++! 

         ! Keyword: RESTART
         ! Added: 30/11/2022 (kr366), last modified: 30/11/2022 (kr366)
         ! Description: Restart simulation from restart file
         ELSE IF (WORD .EQ. 'RESTART') THEN
            RESTARTSIMT = .TRUE.
            CALL READA(RESTARTINPF)
            CALL READA(CONTINUEDUMMY)
            IF (CONTINUEDUMMY.EQ."T") CONTINUESIMT=.TRUE.
            

         ! Keyword: REXMD
         ! Added: 13/12/2022 (kr366), last modified: 13/12/2022 (kr366)
         ! Description: Replica exchange simualtion, either Hamiltonian or Temperature
         ELSE IF (WORD .EQ. 'REXMD') THEN
            REXT= .TRUE.
            CALL READI(NREPLICA)
            CALL READI(NREXSTEPS)
            CALL READA(REXMODE)
            CALL READF(LOWR)
            CALL READF(HIGHR)

         ELSE IF (WORD .EQ. 'RMSD') THEN
            RMSDT = .TRUE.
            CALL READI(NDUMPR)
            IF (NITEMS.GT.2) THEN
               CALL READA(CONTINUEDUMMY)
               IF (CONTINUEDUMMY.EQ."T") ALIGNCONFT = .TRUE.
            END IF

         !+++++++++++++++!   
         ! LETTER S      !
         !+++++++++++++++!

         ! Keyword: SCALEDAT
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Time steps to be used
         ELSE IF (WORD .EQ. 'SCALEDAT') THEN
            CALL READA(SCALEDATNAME)

         !+++++++++++++++!   
         ! LETTER T      !
         !+++++++++++++++!

         ! Keyword: TEMPERATURE
         ! Added: 02/11/2022 (kr366), last modified: 02/11/2022 (kr366)
         ! Description: Temperature to be used
         ELSE IF (WORD .EQ. 'TEMPERATURE') THEN
            CALL READF(TEMP)

         ! Keyword: THERMALISATION
         ! Added: 20/11/2022 (kr366), last modified: 30/11/2022 (kr366)
         ! Description: Thermalisation
         ELSE IF (WORD .EQ. 'THERMALISATION') THEN
            THERMINIT = .TRUE.
            CALL READI(NTHERMALISE)
            CALL READI(NEQUIL)
            IF (NITEMS.GT.3) THEN
               CALL READF(TINIT)
               IF (NITEMS.GT.4) THEN
                  CALL READF(TFINAL)
               ELSE
                  ! set TFINAL to be negative, we then set it to TEMP later
                  TFINAL = -1.0D0
               END IF
            ELSE
               TINIT = 0.0D0
               TFINAL = -1.0D0
            END IF
            

         ! Keyword: THERMALISE_OPTIONS
         ! Added: 20/11/2022 (kr366), last modified: 30/11/2022 (kr366)
         ! Description: Thermalisation settings
         ELSE IF (WORD .EQ. 'THERMALISE_OPTIONS') THEN
            CALL READI(NCENTRE)
            CALL READI(NRMANG)
            CALL READI(NRESCALE)  
            CALL READI(NEQDUMPE)        

         ! Keyword: TIMESTEP
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Time steps to be used
         ELSE IF (WORD .EQ. 'TIMESTEP') THEN
            CALL READF(DT)

         ! Keyword: TOPOLOGY
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Topology to be used
         ELSE IF (WORD .EQ. 'TOPOLOGY') THEN
            CALL READA(TOPNAME)

         !+++++++++++++++!   
         ! LETTER U      !
         !+++++++++++++++!

         !+++++++++++++++!   
         ! LETTER V      !
         !+++++++++++++++!     
         
         !+++++++++++++++!   
         ! LETTER W      !
         !+++++++++++++++!     
         
         !+++++++++++++++!   
         ! LETTER X      !
         !+++++++++++++++!      
               
         !+++++++++++++++!   
         ! LETTER Y      !
         !+++++++++++++++!
         
         !+++++++++++++++!   
         ! LETTER Z      !
         !+++++++++++++++!     
                  
         ELSE
            CALL REPORT('Unrecognized command '//WORD,.TRUE.)
            STOP
         ENDIF
         ! CALL FLUSH(MYUNIT)  
            
      END SUBROUTINE SETKEYS

END MODULE MD_SETUP
