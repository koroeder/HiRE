      SUBROUTINE IO1
      USE commons
      USE genrigid
      USE porfuncs
      IMPLICIT NONE

      INTEGER J1, JP
      IF (DEBUG) THEN
         WRITE(MYUNIT,20) 
20       FORMAT(' io1> Initial coordinates:')
         IF (MPIT) THEN
            WRITE(MYUNIT,30) (COORDS(J1,MYNODE+1),J1=1,3*NATOMS)
         ELSE 
            DO JP=1,NPAR
               WRITE(MYUNIT,30) (COORDS(J1,JP),J1=1,3*NATOMS)
30             FORMAT(3F20.10)
            ENDDO
         ENDIF
      ENDIF
      WRITE(MYUNIT,'(A,I6,A)') ' io1> ', NATOMS,' HiRE beads'
      IF (RIGIDINIT) THEN
         WRITE(MYUNIT, '(I6,A)') NRIGIDBODY,' RIGID BODIES'
         WRITE(MYUNIT, '(I6,A)') (DEGFREEDOMS-6*NRIGIDBODY)/3, ' SINGLE ATOMS'
      ENDIF
      IF (SPARSET) THEN
         WRITE(MYUNIT, '(A,G20.10)') 'io1> Using sparse matrix approach for Hessian, cutoff=',ZERO_THRESH
      ENDIF

      IF (DEBUG) WRITE(MYUNIT,'(A,I6,A)') 'io1> checking the energy of the saved coordinates in the chain'
      IF (FREEZE) THEN
         WRITE(MYUNIT,'(A,I6,A)') 'io1> ', NFREEZE,' atoms will be frozen:'
         DO J1=1,NATOMS
            IF (FROZEN(J1)) WRITE(MYUNIT,'(I6)') J1
         ENDDO
      ENDIF
      IF (HARMONICF) THEN
         WRITE(MYUNIT,'(A,F12.4)') 'io1> harmonically constrained atoms: strength = ', HARMONICSTR
         DO J1=1,NATOMS
            IF (HARMONICFLIST(J1)) WRITE(MYUNIT,'(I6)') J1
         ENDDO
      ENDIF

      IF (RADIUS.EQ.0.0D0) THEN
         RADIUS=2.0D0+(3.0D0*NATOMS/17.77153175D0)**(1.0D0/3.0D0)
         RADIUS=RADIUS*10.0D0
      ENDIF
      IF (NPAR.GT.1) THEN
         WRITE(MYUNIT,'(I2,A)') NPAR,' parallel runs'
      ENDIF
      WRITE(MYUNIT,'(A,G20.10)') 'Sloppy quench tolerance for RMS gradient ',BQMAX
      DO JP=1,NPAR
         IF (FIXBOTH(JP)) THEN
            WRITE(MYUNIT,'(A,I3,A,F12.4,A,2F12.4,A)') &
     &                 'In run ',JP,' temperature=',TEMP(JP),' step size and angular threshold=', &
     &                  STEP(JP),ASTEP(JP),' all fixed'
         ELSE IF (FIXSTEP(JP)) THEN
            WRITE(MYUNIT,'(A,I3,A,2F12.4)') 'In run ',JP,' step size and angular threshold fixed at ', &
     &                                    STEP(JP),ASTEP(JP)
            IF (.NOT.FIXTEMP(JP)) THEN
               WRITE(MYUNIT,'(A,F12.4,A,F12.4)') &
     &                    'Temperature will be adjusted for acceptance ratio ',ACCRAT(JP),' initial value=',TEMP(JP)
            ELSE
               WRITE(MYUNIT,'(A,I1,A,G12.4)') 'In run ',JP,' temperature will be fixed - see below for value'
            ENDIF
         ELSE 
            WRITE(MYUNIT,'(A,I3,A,G12.4)') 'In run ',JP,' temperature will be fixed - see below for value'
            WRITE(MYUNIT,'(A,F12.4,A,2F12.4)') 'Step size and angular threshold will be adjusted for acceptance ratio ', &
     &                ACCRAT(JP),' initial values=',STEP(JP),ASTEP(JP)
         ENDIF
      ENDDO 
      
      WRITE(MYUNIT,'(A)') 'Configuration will be reset to quench geometry'
      WRITE(MYUNIT,'(A)') 'Sampling using Boltzmann weights'      

      WRITE (MYUNIT,'(A)') 'Nocedal LBFGS minimisation'
      WRITE(MYUNIT,'(A,I6)') 'Number of updates before reset in LBFGS=',MUPDATE
      WRITE(MYUNIT,'(A,F20.10)') 'Maximum step size=',MAXBFGS
      WRITE(MYUNIT,'(A,G12.4)') 'Guess for initial diagonal elements in LBFGS=',DGUESS

      WRITE(MYUNIT,'(A,G20.10)') 'Final quench tolerance for RMS gradient ',CQMAX
      WRITE(MYUNIT,'(A,F20.10)') 'Energy difference criterion for minima=',ECONV
      WRITE(MYUNIT,'(A,I5,A,I5)') 'Maximum number of iterations: sloppy quenches ',MAXIT,' final quenches ',MAXIT2
      WRITE(MYUNIT,120) MCSTEPS, TFAC
120         FORMAT('GMIN Run: ',I9,' steps with temperature scaled by ',E15.8)
      IF (DEBUG) THEN
         WRITE(MYUNIT,160) 
160      FORMAT('Debug printing is on')
      ENDIF

      WRITE(MYUNIT, '(A,G20.10)') 'Maximum allowed energy rise during a minimisation=',MAXERISE

      IF (RESTORET) THEN
         WRITE(MYUNIT,'(A,A)') 'Restoring GMIN run from file ',TRIM(ADJUSTL(DUMPFILE))
      ENDIF

      RETURN
      END
