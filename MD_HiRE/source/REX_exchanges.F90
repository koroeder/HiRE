MODULE EXCHANGES
   USE NUMKIND
   USE MD_COMMONS, ONLY: NREPLICA, TASKID
   !> we use alternate exchanges
   !> for odd: 1<->2 3<->4 etc
   !> for even: 2<->3 4<->5 etc
   !> number of rounds of exchnages so far
   INTEGER :: NROUNDEX = 0
   !> number of total successful exchanges
   INTEGER :: NEXCHANGES = 0
   !> current id for node
   INTEGER :: MYCURRENTID
   !> contains the current order of replicas (the first id has the lowest T/highest lambda and the last id the highest T/lowest lambda)
   INTEGER, ALLOCATABLE :: CURRENT_ORDER(:)
   CONTAINS

      SUBROUTINE SELECT_EXCHANGES()
#ifdef MPI
         USE MD_COMMONS, ONLY:  REXMODE, NTASKS, TASKID, MYUNIT
         INTEGER :: NPAIRS !number of replica pairs
         INTEGER :: FIRSTREP !id of first replica
         LOGICAL :: ACTIVEREPT(NREPLICA) !which replicas are active
         LOGICAL :: INITIATE(NREPLICA) !which replicas initiate exchanges
         LOGICAL :: ODDT
         INTEGER :: IDREP, J, NTHISTIME
         LOGICAL :: MEINITIATET, MEACTIVET, EXCHANGEDT
         INCLUDE 'mpif.h'
         INTEGER MPISTATUS(MPI_STATUS_SIZE)         
         INTEGER :: ERR_CODE_MPI, REMD_TAG

         ! increase the round number and set whether this is an odd or even round
         NROUNDEX = NROUNDEX + 1
         ODDT = .FALSE.
         IF (MOD(NROUNDEX,2).EQ.1) ODDT=.TRUE.
         ! set active and initiators to false for all
         ACTIVEREPT(1:NREPLICA) = .FALSE.
         INITIATE(1:NREPLICA) = .FALSE.
         MEACTIVET = .FALSE.
         MEINITIATET = .FALSE.

         IF (TASKID.EQ.0) THEN
            !identify the number of pairs and the first rep
            IF (MOD(NREPLICA,2).EQ.0) THEN
               IF (ODDT) THEN
                  NPAIRS = NREPLICA/2
               ELSE
                  NPAIRS = (NREPLICA-2)/2
               END IF
            ELSE
               NPAIRS = (NREPLICA-1)/2
            END IF
            IF (ODDT) THEN
               FIRSTREP = 1
            ELSE
               FIRSTREP = 2
            END IF
            DO J=FIRSTREP,2*NPAIRS+FIRSTREP-1
               ! which core is the Jth replica on
               IDREP = CURRENT_ORDER(J) + 1
               ! set this one to active
               ACTIVEREPT(IDREP) = .TRUE.
               ! check whether this is an initiator (lower replica)
               IF (ODDT) THEN
                  IF (MOD(J,2).EQ.1) THEN
                     INITIATE(IDREP) = .TRUE.
                  END IF
               ELSE
                  IF (MOD(J,2).EQ.0) THEN
                     INITIATE(IDREP) = .TRUE.
                  END IF
               END IF
            ENDDO
         END IF
         WRITE(MYUNIT,'(A)') " rexmd> Broadcasting the details about exchanges to be tried" 
         !now broad cast the data to all nodes
         CALL MPI_BCAST(ACTIVEREPT,NREPLICA,MPI_LOGICAL,0,MPI_COMM_WORLD,ERR_CODE_MPI)
         CALL MPI_BCAST(INITIATE,NREPLICA,MPI_LOGICAL,0,MPI_COMM_WORLD,ERR_CODE_MPI)

         ! at this stage every node needs to check two things:
         ! a) Am I active? -> MEACTIVET
         ! b) If so, do I initiate? -> MEINITIATET
         IF (ACTIVEREPT(TASKID+1)) THEN
            MEACTIVET = .TRUE.
            IF (INITIATE(TASKID+1)) THEN
               MEINITIATET = .TRUE.
            END IF
         END IF
         ! we can now exclude all replicas that are not active - we will have an MPI_Barrier after this so they won't run away
         IF (MEACTIVET) THEN
            IF (REXMODE.EQ.'T') THEN
               WRITE(MYUNIT, '(A)') " rexmd> Attempting T-REX exchange"
               CALL PASS_DATA_FOR_EXCHANGE_T(MEINITIATET,EXCHANGEDT)
            ELSE IF (REXMODE.EQ.'H') THEN
               WRITE(*,*) "H_REX not implemented"
               STOP
               !CALL PASS_DATA_FOR_EXCHANGE_H(MEINITIATET,EXCHANGEDT)
            END IF
         END IF

         CALL UPDATE_CURR_ORDER()
         IF (TASKID.EQ.0) THEN
            NTHISTIME = 0
            IF (EXCHANGEDT) NTHISTIME = NTHISTIME + 1
            DO J=2,NREPLICA
               CALL MPI_RECV(EXCHANGEDT,1,MPI_LOGICAL,J-1,REMD_TAG,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
               IF (EXCHANGEDT) NTHISTIME = NTHISTIME + 1
            ENDDO
            NEXCHANGES = NEXCHANGES + NTHISTIME/2
            WRITE(*,*) " sel_exchanges> ", NTHISTIME/2, " exchanges this step, in total ", NEXCHANGES, " up to now"
         ELSE
            CALL MPI_SEND(EXCHANGEDT,1,MPI_LOGICAL,0,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
         END IF
         ! TODO: add arguments to MPI barrier - likely need to add MPI COMM RANK calls in thsi routine and in PASS_DATA ...
         ! TODO: write subroutine to H-REX data passing
         ! TODO: add pre-compiler statements and change logging and setup in init routines
         ! TODO: add current_order container to init - should be just the taskids to start
         ! TODO: need to add additional definitions for MPI variables
         ! TODO: add debug printing to test routines
         ! TODO: add exchange steps into MD simulation workflow
         ! TODO: add termination for MD and error handling for replicas into simulation (can we gracefully terminate?)
         ! TODO: in init, assign MYCURRENTID from this module
         CALL MPI_BARRIER(MPI_COMM_WORLD,ERR_CODE_MPI)
#endif
      END SUBROUTINE SELECT_EXCHANGES

      SUBROUTINE PASS_DATA_FOR_EXCHANGE_H(MEINITIATET,EXCHANGEDT)
         USE MD_COMMONS, ONLY: EPOT, TEMP, VEL, NTASKS, COORDS, LAMBDA, NATOMS
         USE HIRE_INTERFACE, ONLY: HIRE_ENERGY_GRAD, SET_UNIV_SCALING
         LOGICAL, INTENT(IN) :: MEINITIATET
         LOGICAL, INTENT(OUT) :: EXCHANGEDT
#ifdef MPI
         INCLUDE 'mpif.h'
         INTEGER MPISTATUS(MPI_STATUS_SIZE)
         REAL(KIND=REAL64) :: U11, U12, U21, U22, T1, T2, L1, L2, PROB, DUMMY
         REAL(KIND=REAL64) :: X1(3*NATOMS), X2(3*NATOMS)
         INTEGER :: OTHERREP, J
         INTEGER :: ERR_CODE_MPI, REMD_TAG
         LOGICAL :: SWITCHT

         ! need to get lambda (L1,L2) from both and coords (X1,X2) from both
         ! set scaling to L1 and get energies for X1 and X2 (U11 and U12)
         ! set scaling to L2 and get energies for X1 and X2 (U21 and U22)

         ! calculate the probability and apply the acceptance/rejection criterion

         ! if accepted:
         ! exchange lambdas and update the order as for T-REX

         ! if rejected:
         ! make sure to reset scaling

#endif 
      END SUBROUTINE PASS_DATA_FOR_EXCHANGE_H

         

      SUBROUTINE PASS_DATA_FOR_EXCHANGE_T(MEINITIATET,EXCHANGEDT)
         USE MD_COMMONS, ONLY: MYUNIT, EPOT, TEMP, VEL, NTASKS, TINIT, TFINAL
         LOGICAL, INTENT(IN) :: MEINITIATET
         LOGICAL, INTENT(OUT) :: EXCHANGEDT
#ifdef MPI
         INCLUDE 'mpif.h'
         INTEGER MPISTATUS(MPI_STATUS_SIZE)
         REAL(KIND=REAL64) :: U1, U2, T1, T2, PROB, DUMMY, RAND, DPRAND
         INTEGER :: OTHERREP, J, I
         INTEGER :: ERR_CODE_MPI, REMD_TAG
         LOGICAL :: SWITCHT 

         ! the initiator receives the information needed for the acceptance/rejection criterion
         IF (MEINITIATET) THEN
            OTHERREP = CURRENT_ORDER(MYCURRENTID+1)
            U1 = EPOT
            T1 = TEMP
            ! set unique tags for each transfer
            REMD_TAG = (TASKID+1)*1000
            CALL MPI_RECV(U2,1,MPI_DOUBLE,OTHERREP,REMD_TAG,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
            CALL MPI_RECV(T2,1,MPI_DOUBLE,OTHERREP,REMD_TAG+1,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
         ELSE
            OTHERREP = CURRENT_ORDER(MYCURRENTID-1)
            ! set the same tag as the receiver
            REMD_TAG = (OTHERREP+1)*1000
            CALL MPI_SEND(EPOT,1,MPI_DOUBLE,OTHERREP,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
            CALL MPI_SEND(TEMP,1,MPI_DOUBLE,OTHERREP,REMD_TAG+1,MPI_COMM_WORLD,ERR_CODE_MPI)
         END IF
         CALL MPI_BARRIER(MPI_COMM_WORLD,ERR_CODE_MPI)
         WRITE(MYUNIT,'(A)') " passdataT> Sent and received data to determine exchanges" 

         ! the value for the probability is given by:
         ! P(1<->2) = min(1, exp[(1/kT1 - 1/kT2)(U1 - U2)]), where Ui is the potential energy of i

         IF (MEINITIATET) THEN
            DUMMY = (1.0/T1 - 1.0/T2)*(U1 - U2)
            PROB = MIN(1.0,EXP(DUMMY))
            RAND = DPRAND()
            WRITE(*,*) " rex> Exchanging ", TASKID, " with ", OTHERREP, &
                       "      DUMMY: ", DUMMY, " ,EXP(DUMMY): ", EXP(DUMMY), &
                       "      Porb: ", PROB, "and random number: ", RAND
            ! accept exchange
            IF (RAND.LT.PROB) THEN
               SWITCHT = .TRUE.
            ! or reject it
            ELSE
               SWITCHT = .FALSE.
            END IF
            ! then pass the status to the other replica (the tags are still set
            ! and should work as there is an MPI_Barrier before this section)
            CALL MPI_SEND(SWITCHT,1,MPI_LOGICAL,OTHERREP,REMD_TAG,MPI_COMM_WORLD,ERR_CODE_MPI)
         ELSE
            CALL MPI_RECV(SWITCHT,1,MPI_LOGICAL,OTHERREP,REMD_TAG,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
         END IF

         ! if we exchange the replicas, we switch their temperature, the order in the current list, and then rescale
         IF (SWITCHT) THEN
            WRITE(MYUNIT,'(A)') " passdataT> exchanging replicas"
            EXCHANGEDT = .TRUE.
            TINIT = TEMP
            CALL MPI_SEND(TEMP,1,MPI_DOUBLE,OTHERREP,REMD_TAG+1,MPI_COMM_WORLD,ERR_CODE_MPI)
            CALL MPI_RECV(TFINAL,1,MPI_DOUBLE,OTHERREP,REMD_TAG+1,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
            TEMP = TFINAL
            !scale the velocities
            VEL(:) = SQRT(TFINAL/TINIT)*VEL(:)
            IF (MEINITIATET) THEN
               WRITE(*,*) " T_REX> exchanged replicas ", MYCURRENTID, " and ", MYCURRENTID+1 , " with temperatures ", TINIT, " and ", TFINAL
               MYCURRENTID = MYCURRENTID + 1
            ELSE 
               MYCURRENTID = MYCURRENTID - 1
            END IF
         ELSE
            WRITE(MYUNIT,'(A)') " passdataT> not exchanging replicas"
            EXCHANGEDT = .FALSE.
         END IF
         WRITE(MYUNIT,'(A)') " "
#endif
      END SUBROUTINE PASS_DATA_FOR_EXCHANGE_T


      SUBROUTINE UPDATE_CURR_ORDER()
         USE MD_COMMONS, ONLY: NREPLICA, TASKID, NTASKS
         IMPLICIT NONE
#ifdef MPI
         INTEGER :: NEWPOS, J
         INTEGER :: NEWORDER(NREPLICA)
         INTEGER :: ERR_CODE_MPI
         INCLUDE 'mpif.h'
         INTEGER MPISTATUS(MPI_STATUS_SIZE)

         IF (TASKID.EQ.0) THEN
            NEWORDER(1:NREPLICA) = -1
            NEWORDER(MYCURRENTID) = 0
            DO J=1,NTASKS-1
               CALL MPI_RECV(NEWPOS,1,MPI_INTEGER,J,J,MPI_COMM_WORLD,MPISTATUS,ERR_CODE_MPI)
               NEWORDER(NEWPOS) = J 
            END DO
            CURRENT_ORDER = NEWORDER
         ELSE
            CALL MPI_SEND(MYCURRENTID,1,MPI_INTEGER,0,TASKID,MPI_COMM_WORLD,ERR_CODE_MPI)
         END IF
         CALL MPI_BCAST(CURRENT_ORDER,NREPLICA,MPI_INTEGER,0,MPI_COMM_WORLD,ERR_CODE_MPI)
#endif
      END SUBROUTINE UPDATE_CURR_ORDER
! get dedicated MPI function for everything - this will make things easier

! key question: keep process id with same replica or keep T in order

! will keep process id with same replica -> book keeping needed (order o T for all replicas)
! allows to keep log files associated with each individual node
! on eschange - only need to exchange Ts and kambdas

! what needs to be exhanged:
! 1. set odd or even exchange step (can use an overall counter rather than logical ...)
! 2. decide on all pairs to be tested for exchange
! 3. send all replicas that might be exchanged a signal
! 4. now send all involved replicas a signal if they are the initiator (lower T/higher lambda)
! 5. a) initiator receive:
!      - T2 and U2 for T-REX      
!      - T2, x2 and lambda for H-REX
!    b) initiator check acceptance for T-REX
!    c) for H-REX: calculate U1 and U2 with each x1 and x2, then calculate acceptance
! 6. initiator send result to head
! 7. head call exchange of replicas

! exchange:
! for T-REX
! make higher and lower T available to both replicas
! exchange velocities and coordinates
! scale velocities to new T (see scaling for initiation, sqrt of Tnew/Told)

! for H-REX
! make lambda available to both
! exchange velocities and coordinates
! is there a need to rescle velocities? Or will they be auto adjusted? Looks like no adjustment is done here inn the old OPEP/HiRE simulateur


!! Acceptance/rejection
!
!T-REX: P(1<->2) = min(1, exp[(1/kT1 - 1/kT2)(U1 - U2)]), where Ui is the potential energy of i
!
!H-REX: P(1<->2) = min(1, exp[(1/kT1 - 1/kT2){(U1(x2)-U1(x1)) + (U2(x1)-U2(x2))}])
!


END MODULE EXCHANGES
