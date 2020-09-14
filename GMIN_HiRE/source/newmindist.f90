!   Finds the minimum distance between two geometries.
!   Geometry in RA should not change. RB is returned as the
!   closest geometry to RA if PRESERVET is .FALSE.
!
!   New analytic method based on quaterions from
!   Kearsley, Acta Cryst. A, 45, 208-210, 1989.
!
!   Modified for general angle-axis 30/01/12
!  

SUBROUTINE NEWMINDIST(RA,RB,NATOMS,DIST,DEBUG,RMAT)
USE COMMONS,ONLY : MYUNIT
USE PREC
IMPLICIT NONE
INTEGER J1, NATOMS, NSIZE, INFO, JMIN
INTEGER,PARAMETER :: LWORK=12
REAL(KIND = REAL64) RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), XM, YM, ZM, XP, YP, ZP, &
  &              DIAG(4), TEMPA(LWORK), RMAT(3,3), MINV, Q1, Q2, Q3, Q4, CMXA, CMYA, CMZA, CMXB, CMYB, CMZB

REAL(KIND = REAL64), ALLOCATABLE :: XA(:), XB(:)
LOGICAL DEBUG

REAL(KIND = REAL64) XSHIFT, YSHIFT, ZSHIFT

   ALLOCATE(XA(3*NATOMS),XB(3*NATOMS))
   NSIZE=NATOMS
   XA(1:3*NATOMS)=RA(1:3*NATOMS)
   XB(1:3*NATOMS)=RB(1:3*NATOMS)

   ! Move centre of coordinates of XA and XB to the origin.
   DO J1=1,NSIZE
      CMXA=CMXA+XA(3*(J1-1)+1)
      CMYA=CMYA+XA(3*(J1-1)+2)
      CMZA=CMZA+XA(3*(J1-1)+3)
   ENDDO
   CMXA=CMXA/NSIZE; CMYA=CMYA/NSIZE; CMZA=CMZA/NSIZE
   DO J1=1,NSIZE
      XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-CMXA
      XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-CMYA
      XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-CMZA
   ENDDO
   DO J1=1,NSIZE
      CMXB=CMXB+XB(3*(J1-1)+1)
      CMYB=CMYB+XB(3*(J1-1)+2)
      CMZB=CMZB+XB(3*(J1-1)+3)
   ENDDO
   CMXB=CMXB/NSIZE; CMYB=CMYB/NSIZE; CMZB=CMZB/NSIZE
   DO J1=1,NSIZE
      XB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB
      XB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB
      XB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB
   ENDDO

   XSHIFT=0.0D0; YSHIFT=0.0D0; ZSHIFT=0.0D0

!  The formula below is not invariant to overall translation because XP, YP, ZP
!  involve a sum of coordinates! We need to have XA and XB coordinate centres both
!  at the origin!!
   QMAT(1:4,1:4)=0.0D0
!  PRINT *,'XA:'
!  PRINT '(6G20.10)',XA(1:3*NATOMS)
!  PRINT *,'XB:'
!  PRINT '(6G20.10)',XB(1:3*NATOMS)
   DO J1=1,NSIZE
      XM=XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
      YM=XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
      ZM=XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
      XP=XA(3*(J1-1)+1)+XB(3*(J1-1)+1)
      YP=XA(3*(J1-1)+2)+XB(3*(J1-1)+2)
      ZP=XA(3*(J1-1)+3)+XB(3*(J1-1)+3)
      QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
      QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
      QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
      QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
      QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
      QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
      QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
      QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
      QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
      QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
   ENDDO
   QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3); QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)

   CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,LWORK,INFO)
   IF (INFO.NE.0) WRITE(MYUNIT,'(A,I6,A)') 'newmindist> WARNING - INFO=',INFO,' in DSYEV'

   MINV=1.0D100
   DO J1=1,4
      IF (DIAG(J1).LT.MINV) THEN
         JMIN=J1
         MINV=DIAG(J1)
      ENDIF
   ENDDO
   IF (MINV.LT.0.0D0) THEN
      IF (ABS(MINV).LT.1.0D-6) THEN
         MINV=0.0D0
      ELSE
         WRITE(MYUNIT,'(A,G20.10,A)') 'newmindist> WARNING MINV is ',MINV,' change to absolute value'
         MINV=-MINV
      ENDIF
   ENDIF
   DIST=SQRT(MINV)

   IF (DEBUG) WRITE(MYUNIT,'(A,G20.10,A,I6)') 'newmindist> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
   Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
!
! RMAT will contain the matrix that maps RB onto the best correspondence with RA
!
   RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
   RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
   RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
   RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
   RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
   RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
   RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
   RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
   RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2


   ! Translate the RB coordinates to the centre of coordinates of RA. 
   DO J1=1,NATOMS
      RB(3*(J1-1)+1)=RB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
      RB(3*(J1-1)+2)=RB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
      RB(3*(J1-1)+3)=RB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
   ENDDO
   CALL NEWROTGEOM(NSIZE,RB,RMAT,CMXA,CMYA,CMZA)

DEALLOCATE(XA,XB)

END SUBROUTINE NEWMINDIST


SUBROUTINE NEWROTGEOM(NATOMS,COORDS,ROTMAT,CX,CY,CZ)
USE PREC
IMPLICIT NONE
INTEGER I, J, K, NATOMS
REAL(KIND = REAL64) COORDS(*), R1, R0(3), ROTMAT(3,3), CX, CY, CZ

DO I=1,NATOMS
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+ROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE NEWROTGEOM


