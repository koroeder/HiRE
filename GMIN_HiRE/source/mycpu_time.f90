
SUBROUTINE MYCPU_TIME(POO)
USE COMMONS
IMPLICIT NONE
REAL(KIND = REAL64) MYTIME, POO

CALL CPU_TIME(MYTIME)
POO=MYTIME ! without this extra assignment NAG f95 returns a random number!
           ! saving TSTART is necessary for PG compiler, where initial time is
           ! not zero!
RETURN
END SUBROUTINE MYCPU_TIME
