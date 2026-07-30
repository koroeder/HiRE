PROGRAM GMIN
!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2019 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
      USE NOA     
      USE COMMONS
      USE PORFUNCS
      USE GENRIGID, only : RIGIDINIT, GENRIGID_READ_FROM_FILE, DEGFREEDOMS

      IMPLICIT NONE
#ifdef MPI
      INCLUDE 'mpif.h'
      INTEGER :: MPIERR
#endif
      CHARACTER(LEN=130) :: ISTR
      CHARACTER(LEN=10) :: JSTR
      CHARACTER(LEN=256) :: MYFILENAME

      INTEGER J1, JP, GETUNIT
      LOGICAL LOPEN

      CALL CPU_TIME(TSTART)
      CALL READ_CMD_ARGS
#ifdef MPI
      CALL MPI_INIT(MPIERR)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NPAR,MPIERR)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYNODE,MPIERR)
      MYUNIT=22980+MYNODE
      WRITE (ISTR, '(I10)') MYNODE+1
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') INFILE='output'
      IF (NPAR.GT.1) THEN 
         MYFILENAME=TRIM(ADJUSTL(INFILE))//"."//TRIM(ADJUSTL(ISTR))
      ELSE
         MYFILENAME=TRIM(ADJUSTL(INFILE))
      ENDIF
      OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      WRITE(MYUNIT, '(A,I10,A,I10)') "Starting parallel execution: Processor", mynode+1, " of ",NPAR
#else
      NPAR=1
      MYNODE=0
      IF (TRIM(ADJUSTL(INFILE)).EQ.'') THEN
         MYUNIT=6
      ELSE
         MYUNIT=22979+1
         MYFILENAME=TRIM(ADJUSTL(INFILE))
         OPEN(MYUNIT,FILE=MYFILENAME, STATUS="unknown", form="formatted")
      ENDIF
      WRITE(MYUNIT,'(A)') 'Starting serial execution'
#endif

      CALL COUNTATOMS()
      NATOMSALLOC=NUMBER_OF_ATOMS
      NATOMS=NUMBER_OF_ATOMS

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> A ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF
      CALL KEYWORD

      IF (RIGIDINIT) THEN
          CALL GENRIGID_READ_FROM_FILE ()
          WRITE(MYUNIT, '(A,I15)') " genrigid> rbodyconfig used to specifiy rigid bodies, degrees of freedom now ", DEGFREEDOMS
      END IF

      IF (FEBHT) THEN
         IF (.NOT.(ALLOCATED(ATMASS))) THEN
            ALLOCATE(ATMASS(NATOMSALLOC))
            ATMASS(:)=1.D0            
         ENDIF 
         SYMFCTR=LOG(2.D0) ! this is the factor of 2 for inversion
         DO J1=1, NATOMSALLOC
            SYMFCTR=SYMFCTR+LOG(1.0D0*J1)
         ENDDO
      ENDIF


      IF (.NOT.ALLOCATED(ATMASS)) THEN
         WRITE(*,'(A)') 'main> No atomic masses allocated'
      ENDIF

      ALLOCATE(LNFAC(NATOMSALLOC)) ! store ln factorial values
      LNFAC(1)=0.0D0
      DO J1=2,NATOMS
         LNFAC(J1)=LNFAC(J1-1)+LOG(1.0D0*J1)
      ENDDO

      INQUIRE(UNIT=1,OPENED=LOPEN)
      IF (LOPEN) THEN
         WRITE(*,'(A,I2,A)') 'main> B ERROR *** Unit ', 1, ' is not free '
         STOP
      ENDIF


      ALLOCATE(FF(NSAVE),QMIN(MAX(NSAVE,1)))
      ALLOCATE(QMINP(NSAVE,3*NATOMSALLOC))
      ALLOCATE(QMINNATOMS(NSAVE))
      ALLOCATE(QMINT(NSAVE,NATOMSALLOC))
      ALLOCATE(NPCALL_QMIN(NSAVE))

      QMINP(1:NSAVE,1:3*NATOMSALLOC)=0.0D0 ! to prevent reading from uninitialised memory
      QMINT(1:NSAVE,1:NATOMSALLOC)=1 ! to prevent reading from uninitialised memory
      QMINNATOMS(1:NSAVE)=NATOMS ! to prevent reading from uninitialised memory
      COORDSO(1:3*NATOMSALLOC,1:NPAR)=0.0D0 ! to prevent reading from uninitialised memory
      FF(1:NSAVE)=0 ! to prevent reading from uninitialised memorY


      IF (DUMP_MARKOV) THEN
         ! dump.1.xyz is partly filled with control characters for no apparent reason.
         ! Suspect a compiler or mpi bug? DJW
         ! dump.1.xyz is fine with debug compilation, so it is a compiler bug!
         ALLOCATE(DUMPXYZUNIT(NPAR))
         DO J1=1,NPAR
            WRITE (JSTR,'(I6)') J1
            ISTR='dump.' // TRIM(ADJUSTL(JSTR)) // '.xyz'
            DUMPXYZUNIT(J1)=GETUNIT()
            OPEN(UNIT=DUMPXYZUNIT(J1),FILE=TRIM(ADJUSTL(ISTR)),STATUS='UNKNOWN')
         ENDDO
      ENDIF


!     TRACKDATA keyword prints the energy and markov energy to files for viewing during a run.
     IF (TRACKDATAT) THEN
         MYEUNIT=4000+MYNODE
         MYMUNIT=6000+MYNODE
         MYBUNIT=8000+MYNODE
         IF (FEBHT) FE_FILE_UNIT=12000+MYNODE
         IF (SAXST) THEN 
            MYSEUNIT=4500+MYNODE
            MYSMUNIT=6500+MYNODE
         ENDIF
         IF (NPAR.GT.1) THEN
            IF (RESTORET) THEN
               OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               IF (FEBHT) THEN
                  OPEN(FE_FILE_UNIT,FILE="free_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               ENDIF
               IF (SAXST) THEN
                  OPEN(MYSEUNIT,FILE="saxs_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
                  OPEN(MYSMUNIT,FILE="saxs_markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               ENDIF 
            ELSE
               OPEN(MYEUNIT,FILE="energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYMUNIT,FILE="markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYBUNIT,FILE="best."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               IF (FEBHT) THEN
                  OPEN(FE_FILE_UNIT,FILE="free_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               ENDIF
               IF (SAXST) THEN
                  OPEN(MYSEUNIT,FILE="saxs_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
                  OPEN(MYSMUNIT,FILE="saxs_markov."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               ENDIF 
            ENDIF
         ELSE
            IF (RESTORET) THEN
               OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               IF (FEBHT) THEN
                  OPEN(FE_FILE_UNIT,FILE="free_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               ENDIF
               IF (SAXST) THEN
                  OPEN(MYSEUNIT,FILE="saxs_energy",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
                  OPEN(MYSMUNIT,FILE="saxs_markov",STATUS='UNKNOWN',FORM='FORMATTED',POSITION='APPEND')
               ENDIF
            ELSE
               OPEN(MYEUNIT,FILE='energy',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYMUNIT,FILE='markov',STATUS='UNKNOWN',FORM='FORMATTED')
               OPEN(MYBUNIT,FILE='best',STATUS='UNKNOWN',FORM='FORMATTED')
               IF (FEBHT) THEN
                  OPEN(FE_FILE_UNIT,FILE="free_energy."//trim(adjustl(istr)),STATUS='UNKNOWN',FORM='FORMATTED')
               ENDIF
               IF (SAXST) THEN
                  OPEN(MYSEUNIT,FILE="saxs_energy",STATUS='UNKNOWN',FORM='FORMATTED')
                  OPEN(MYSMUNIT,FILE="saxs_markov",STATUS='UNKNOWN',FORM='FORMATTED')
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      CALL FLUSH(6)
      CALL IO1


      IF (CENT) THEN
         DO J1=1,NPAR
            CALL CENTRE2(COORDS(1:3*NATOMS,J1))
         ENDDO
      ENDIF

      DO JP=1,NPAR
         NQ(JP)=1
      ENDDO
      DO J1=1,NSAVE
         QMIN(J1)=1.0D10
         NPCALL_QMIN(J1)=0
      ENDDO

      CALL INITIALIZATIONS()

      CALL MCRUNS
   
      CALL FLUSH(MYUNIT)
      CLOSE(MYUNIT)
    
! Close file storing free energies
      IF (FEBHT) THEN
         CLOSE(FE_FILE_UNIT)
      END IF
      IF (TRACKDATAT) THEN
         CLOSE(MYEUNIT)
         CLOSE(MYMUNIT)
         CLOSE(MYBUNIT)
      ENDIF

      IF (ALLOCATED(FF)) DEALLOCATE(FF)
      IF (ALLOCATED(QMIN)) DEALLOCATE(QMIN)
      IF (ALLOCATED(QMINP)) DEALLOCATE(QMINP)
      IF (ALLOCATED(QMINT)) DEALLOCATE(QMINT)
      IF (ALLOCATED(QMINNATOMS)) DEALLOCATE(QMINNATOMS)
      IF (ALLOCATED(FIXSTEP)) DEALLOCATE(FIXSTEP)
      IF (ALLOCATED(FIXTEMP)) DEALLOCATE(FIXTEMP)
      IF (ALLOCATED(FIXBOTH)) DEALLOCATE(FIXBOTH)
      IF (ALLOCATED(TEMP)) DEALLOCATE(TEMP)
      IF (ALLOCATED(ACCRAT)) DEALLOCATE(ACCRAT)
      IF (ALLOCATED(STEP)) DEALLOCATE(STEP)
      IF (ALLOCATED(ASTEP)) DEALLOCATE(ASTEP)
      IF (ALLOCATED(OSTEP)) DEALLOCATE(OSTEP)
      IF (ALLOCATED(NQ)) DEALLOCATE(NQ)
      IF (ALLOCATED(EPREV)) DEALLOCATE(EPREV)          
      IF (ALLOCATED(COORDS)) DEALLOCATE(COORDS)
      IF (ALLOCATED(COORDSO)) DEALLOCATE(COORDSO)
      IF (ALLOCATED(FROZEN)) DEALLOCATE(FROZEN)
      IF (ALLOCATED(FROZENRES)) DEALLOCATE(FROZENRES)
      IF (ALLOCATED(HARMONICFLIST)) DEALLOCATE(HARMONICFLIST)
      IF (ALLOCATED(HARMONICR0)) DEALLOCATE(HARMONICR0)

#ifdef MPI
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERR)
       CALL MPI_FINALIZE(MPIERR)
#endif
      STOP 'main> GMIN has terminated normally.'
      END
