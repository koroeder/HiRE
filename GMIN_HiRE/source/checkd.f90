! Testing routine to check the potential, gradient and Hessian calculations

     SUBROUTINE CHECKD(X)

      USE COMMONS
      USE GENRIGID, ONLY: RIGIDINIT, ATOMRIGIDCOORDT, DEGFREEDOMS, TRANSFORMCTORIGID
      USE PREC
      USE MODHESS
      USE HIRE_INTERFACE, ONLY: TERMINATE_HIRE, DUMP_PDB
      IMPLICIT NONE

      INTEGER             :: IVRNO, IVRNO1, IVRNO2, dof, ITERATIONS, MYUNIT2, J2
      REAL(KIND = REAL64) :: X(3*NATOMS), G(3*NATOMS), ENERGY, FM, FP, DFA, DFN, TMPCOORDS(3*NATOMS)
      REAL(KIND = REAL64) :: TIME1, TIME2, RMSF
      LOGICAL             :: GTEST, STEST, CFLAG
      REAL(KIND = REAL64), PARAMETER :: ERRLIM = 1.D-04, DELX = 1.0D-4

      ! allow for rigid bodies
      if (rigidinit) then
         dof = degfreedoms
      else
         dof = 3*natoms
      endif

      ! WRITE(MYUNIT,'(A,G20.10)') 'DELX: ', DELX

      CFLAG = .FALSE.
      STEST = .FALSE.

      IF (CHECKDID == 0) THEN
         GTEST = .FALSE.
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST, .FALSE.)
         WRITE(MYUNIT, *) 'Energy  = ', ENERGY
         RMSF=MAX(SQRT(SUM(G(1:3*NATOMS)**2)/(3*NATOMS)), 1.0D-100 )
         WRITE(MYUNIT,'(A,2G15.7)') "RMS force: ", RMSF,RMS
      ELSEIF (CHECKDID == 1) THEN
         ! Checks gradients
         ! check derivatives wrt atomic positions
         DO IVRNO = 1, DOF
            WRITE(MYUNIT, *) IVRNO
            if (rigidinit.and.atomrigidcoordt) then
               call transformctorigid(x, tmpcoords)
               x(1:degfreedoms) = tmpcoords(1:degfreedoms)
               x(degfreedoms+1:3*natoms) = 0.0d0
               atomrigidcoordt = .false.
            endif
            GTEST = .FALSE.
            X(IVRNO) = X(IVRNO) - DELX
            CALL POTENTIAL(X, G, FM, GTEST, STEST, .FALSE.)

            if (rigidinit.and.atomrigidcoordt) then
               call transformctorigid(x, tmpcoords)
               x(1:degfreedoms) = tmpcoords(1:degfreedoms)
               x(degfreedoms+1:3*natoms) = 0.0d0
               atomrigidcoordt = .false.
            endif

            X(IVRNO) = X(IVRNO) + 2.D0*DELX
            CALL POTENTIAL(X, G, FP, GTEST, STEST, .FALSE.)
     
            if (rigidinit.and.atomrigidcoordt) then
               call transformctorigid(x, tmpcoords)
               x(1:degfreedoms) = tmpcoords(1:degfreedoms)
               x(degfreedoms+1:3*natoms) = 0.0d0
               atomrigidcoordt = .false. 
            endif
 
            GTEST = .TRUE.
            X(IVRNO) = X(IVRNO) - DELX
            CALL POTENTIAL(X, G, ENERGY, GTEST, STEST, .FALSE.)
            DFN = (FP - FM) / (2.D0*DELX)
            IF (ABS(DFN).LT.1.0D-10) DFN = 0.D0
            DFA = G(IVRNO)

            WRITE(MYUNIT, *) 'Gradient numerical  = ', DFN
            WRITE(MYUNIT, *) 'Gradient analytical = ', DFA

            IF (ABS(DFN - DFA) > ERRLIM) WRITE(MYUNIT, '(A,I10,3G20.10)') 'WARNING *** ',IVRNO, DFN, DFA, ABS(DFN-DFA)
         ENDDO

      ELSE IF (CHECKDID == 2) THEN

         IF (.NOT. ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))

         DO IVRNO1 = 1, 3*NATOMS
            DO IVRNO2 = 1, 3*NATOMS
               WRITE(MYUNIT,*) IVRNO1, IVRNO2
               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE., .FALSE.)
               FM   = G(IVRNO2)

               X(IVRNO1) = X(IVRNO1) + 2.D0*DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE., .FALSE.)
               FP   = G(IVRNO2)

               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.TRUE., .FALSE.)
               DFN  = (FP - FM) / (2.D0*DELX)
               DFA  = HESS(IVRNO1,IVRNO2)

               WRITE(MYUNIT, *) 'Hessian numerical  = ', DFN
               WRITE(MYUNIT, *) 'Hessian analytical = ', DFA

               IF (ABS(DFN - DFA) > ERRLIM) WRITE(MYUNIT,*) 'Error:', IVRNO1, IVRNO2, DFN, DFA, ABS(DFN-DFA)
            ENDDO
         ENDDO

      ELSE IF (CHECKDID == 3) THEN
            WRITE(MYUNIT,'(A)') 'checkd > Minimise configuration'
            CALL POTENTIAL(X, G, ENERGY, GTEST, STEST, .FALSE.) 
            WRITE(MYUNIT, *) 'checkd > Initial energy: ', ENERGY
            CALL CPU_TIME(TIME1)
            CALL MYLBFGS(3*NATOMS,MUPDATE,X,.FALSE.,QTESTMAX,CFLAG,POTEL,MAXIT,ITERATIONS,.TRUE.)
            CALL CPU_TIME(TIME2)
            IF (CFLAG) THEN 
               WRITE(MYUNIT,'(A,G20.10,A,I5,A,G12.5,A,F11.1)') &
                            ' E=', POTEL,' steps=',ITERATIONS,' RMS=',RMS, &
                            ' t=',TIME2-TIME1
               CALL DUMP_PDB(3*NATOMS, X,'quench.pdb',.TRUE.)
               OPEN(MYUNIT2,FILE='start.quench',STATUS="unknown",form="formatted")
               DO J2=1,NATOMS
                  WRITE(MYUNIT2,'(3F28.20)') X(3*(J2-1)+1),X(3*(J2-1)+2),X(3*(J2-1)+3)
               ENDDO
               CLOSE(MYUNIT2)
            ELSE
               WRITE(MYUNIT,'(A)') 'checkd > Quench not converged.'
            ENDIF

      ENDIF
      CALL TERMINATE_HIRE()
      STOP
      END SUBROUTINE CHECKD
