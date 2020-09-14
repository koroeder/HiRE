      SUBROUTINE GSAVEIT(EREAL,P,NP)
      USE COMMONS
      USE DEFS_MCRUNS
      IMPLICIT NONE


      INTEGER J1, J2, J3, NP
      REAL(KIND = REAL64) EREAL,P(3*NATOMS)


      ! Save the lowest NSAVE distinguishable configurations.
      DO J1=1,NSAVE
         IF (DABS(EREAL-QMIN(J1)).LT.ECONV) THEN
            ! These are probably the same - but just to make sure we save the lowest.
            IF (EREAL.LT.QMIN(J1)) THEN
               QMINNATOMS(J1)=NATOMS
               QMIN(J1)=EREAL
               DO J2=1,3*NATOMS
                  QMINP(J1,J2)=P(J2)
               ENDDO
            ENDIF
            GOTO 10
         ENDIF
         IF (EREAL.LT.QMIN(J1)) THEN

            J2=NSAVE
20          CONTINUE
      
            IF (NSAVE.GT.1) THEN
               QMIN(J2)=QMIN(J2-1)
               QMINNATOMS(J2)=QMINNATOMS(J2-1)
               FF(J2)=FF(J2-1)
               NPCALL_QMIN(J2)=NPCALL_QMIN(J2-1)
               DO J3=1,3*QMINNATOMS(J2-1)
                  QMINP(J2,J3)=QMINP(J2-1,J3)
               ENDDO
               J2=J2-1
               IF (J2.GE.J1+1) GOTO 20
            ENDIF

            QMIN(J1)=EREAL
            QMINNATOMS(J1)=NATOMS
            FF(J1)=NQ(NP)
            NPCALL_QMIN(J1)=NPCALL
            DO J2=1,3*QMINNATOMS(J1)
               QMINP(J1,J2)=P(J2)
            ENDDO
            GOTO 10
         ENDIF
      ENDDO

10    CONTINUE

      RETURN
      END



