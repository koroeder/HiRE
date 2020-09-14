
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Get PFMIN value for minimum J1 at temperature LOCALFETEMP
! We calculate  ln( n! qrot * qvib * exp(-V/kT) / o ) with o the order of the point group
! qvib can be classical or quantum (usefrqs, needs all frequencies from min.frqs file)
! qrot can be polyatomic or linear
! special cases for single atom and diatomic
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE GETPFMIN(LOCALFETEMP,PFMINLOCAL,LNATOMSMIN,LIXMIN,LIYMIN,LIZMIN,LFVIBMIN,LEMIN,LHORDERMIN,LOCALFRQS)  
      USE COMMONS, ONLY : PLANCK, LNFAC, NATOMS, USEFRQS
      USE PREC
      IMPLICIT NONE
      !Arguments
      INTEGER             :: LNATOMSMIN, LHORDERMIN
      REAL(KIND = REAL64) :: LOCALFETEMP, PFMINLOCAL, LIXMIN, LIYMIN, LIZMIN, LFVIBMIN, LEMIN, LOCALFRQS(3*NATOMS)

      !Variables
      INTEGER :: NAT, NUM_ZERO_EVS, NFRQS, J1
      REAL(KIND = REAL64) :: ROTFACPOLY, LNTEMP, LNITDETFAC, LNFRQFAC, LNTFAC
      REAL(KIND = REAL64) :: FRQFAC, ROTTERM, SYMTERM
      REAL(KIND = REAL64) :: VIBTERM, DUMMY, ITDET
      REAL(KIND = REAL64), PARAMETER :: PI=3.141592654D0 

      NAT=LNATOMSMIN
      ROTFACPOLY=LOG(SQRT(8.0D0*PI*LOCALFETEMP**3)/(PLANCK/(2.0D0*PI))**3)
      LNTEMP=LOG(2.0D0*PI*LOCALFETEMP/PLANCK)
!
! Unit conversions for TIP4P. We only need nonlinear polyatomic rotation/vibration.
! Log( sqrt( kJ/mol / (amu * Angstrom^2) ) ) the 1/(2*pi) factor is already included below
! FETEMP and EREAL and mu are in kJ/mol
! PLANCK should be in Js, normal SI
!

     !FRQFAC=PLANCK/(2.0D0*PI)
      FRQFAC=0.310428D0 ! 6.626*10^(-34)*6.022*10^(23)/(10^3*4.184)*Sqrt[10^3*4.184/(6.022*10^(23)*1.661*10^(-27)*10^(-20))]/(2*Pi)
      LNFRQFAC=30.6491  ! log( sqrt(10^3*4.184/(6.022*10^(23)) / (10^(-20)*1.661*10^(-27))) )  - we have already included 1/(2pi)
      LNITDETFAC=-323.142D0 ! log (amu * Angstrom^2)^3 always polyatomic for HiRE - converts to SI
      LNTFAC=-46.4159D0     ! ln( kcal/mol to J )

      NUM_ZERO_EVS=6
      NFRQS=3*NAT-NUM_ZERO_EVS
      IF (USEFRQS) THEN
         VIBTERM=0.0D0
         DO J1=NUM_ZERO_EVS+1,NFRQS+NUM_ZERO_EVS
            DUMMY=-FRQFAC*SQRT(LOCALFRQS(J1))/(2.0D0*LOCALFETEMP)
            VIBTERM=VIBTERM+DUMMY-LOG(1.0D0-EXP(2.0D0*DUMMY))
         ENDDO
      ELSE
         VIBTERM=NFRQS*(LNTEMP  + LNTFAC - LNFRQFAC) - LFVIBMIN/2.0D0
      ENDIF
      ITDET=LIXMIN*LIYMIN*LIZMIN
      ROTTERM=ROTFACPOLY + 1.5D0*LNTFAC+ LOG(ITDET)/2.0D0  + LNITDETFAC/2.0D0
      SYMTERM=LNFAC(NAT)+LOG(1.0D0/LHORDERMIN)
      PFMINLOCAL= -LEMIN/LOCALFETEMP + SYMTERM + VIBTERM + ROTTERM

      END 
