MODULE MD_SETUP

   CONTAINS
      SUBROUTINE START_TRACKING()
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE MD_COMMONS
         IMPLICIT NONE
         CALL FILE_OPEN("md_energy.log",EUNIT,.TRUE.)
         CALL FILE_OPEN("md_coords.xyz",XUNIT,.TRUE.)
      END SUBROUTINE START_TRACKING

      SUBROUTINE SETUP_POTENTIAL()
         USE MD_COMMONS, ONLY: NATOMS, TOPNAME, SCALEDATNAME, COORDSFILE
         USE HIRE_INTERFACE, ONLY: HIRE_INITIALISE, PASS_HIRE_MASSES
         USE MD_UTILS, ONLY: ALLOC_COMMONS
         USE FILE_UTILS, ONLY: FILE_EXIST, FILE_OPEN
         IMPLICIT NONE
         INTEGER :: J, XUNIT
         !first initialise the HiRE interface
         CALL HIRE_INITIALISE(TOPNAME, SCALEDATNAME, NATOMS)

         ! allocate the relevant arrays
         CALL ALLOC_COMMONS()

         ! get coordinates
         IF (FILE_EXIST(COORDSFILE)) THEN
            CALL FILE_OPEN(COORDSFILE,XUNIT,.FALSE.)
            READ(XUNIT, *) (X(J), J=1,3*NATOMS)
            CLOSE(XUNIT)
         ELSE
            WRITE(MYUNIT,*) " setup> Cannot locate input file for coordinates - ", COORDSFILE
            STOP
         END IF

         ! get particle masses
         CALL PASS_HIRE_MASSES(NATOMS,MASSES)
      END SUBROUTINE SETUP_POTENTIAL

      SUBROUTINE READ_SETTINGS(PARAMUNIT)
         USE MD_COMMONS, ONLY: MYUNIT
         USE FILE_UTILS, ONLY: FILE_OPEN
         USE INPUTMOD, ONLY: INPUTKW
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: PARAMUNIT
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
         ! done with reading keywords
         WRITE(MYUNIT, '(A)') " read_keywords> Completed reading simdata"
      END SUBROUTINE READ_SETTINGS

      SUBROUTINE SETKEYS(WORD)
         USE INPUTMOD
         USE MD_COMMONS                  ! global variables
         USE FILE_UTILS, ONLY: FILE_EXIST
        
         IMPLICIT NONE
         INTEGER :: J1
         CHARACTER(25), INTENT(IN) :: WORD
        
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

         ! Keyword: MDSTEPS
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Number of MD steps to be taken
         ELSE IF (WORD .EQ. 'MDSTEPS') THEN
            CALL READI(MDSTEPS) 

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
         ! Description: Time steps to be used
         ELSE IF (WORD .EQ. 'TEMPERATURE') THEN
            CALL READF(TEMP)

         ! Keyword: TIMESTEP
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Time steps to be used
         ELSE IF (WORD .EQ. 'TIMESTEP') THEN
            CALL READF(DT)

         ! Keyword: TOPOLOGY
         ! Added: 01/11/2022 (kr366), last modified: 01/11/2022 (kr366)
         ! Description: Time steps to be used
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
         CALL FLUSH(MYUNIT)  
            
      END SUBROUTINE SETKEYS

END MODULE MD_SETUP