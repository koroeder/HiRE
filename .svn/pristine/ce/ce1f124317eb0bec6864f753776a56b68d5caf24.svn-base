
      SUBROUTINE EIG(A,B,L,N,N1)
!
! DIAGONALIZATION BY THE JACOBI METHOD.
! A - MATRIX TO BE DIAGONALIZED (eigenvalues returned in diagonal
!        elements of A).  If you want to save A, you must do this before
!        calling EIG.  Set N to the same value as L.
! B - EIGENVECTORS
! L - DIMENSION OF A AND B
! N - SIZE OF SUBMATRIX USED
! N1 - A FLAG INDICATING WHETHER THE EIGENVECTORS AND
!      EIGENVALUES ARE TO BE REORDERED.
!

      IMPLICIT NONE
      INTEGER L,N,N1,JJ,IOFF,I,IM1,MU,MM,JI,II,J
      DOUBLE PRECISION A(L,L),B(L,L),W2,W1,C,T,ALP,ALN,D,SUM,S,DIFF,R,Q,P,TOL2,TOL,ZER,ONE
      DATA ZER/0.D00/,ONE/1.D00/

      TOL=1.D-14
      TOL2=1.D-10
      JJ=0
      IOFF=0
      B(1,1)=ONE
      IF(N.EQ.1) RETURN
      DO I=2,N
        IM1=I-1
        DO J=1,IM1
          B(I,J)=ZER
          B(J,I)=ZER
        END DO
        B(I,I)=ONE
      END DO
!
! FIRST SEE IF MATRIX IS ALREADY DIAGONAL- IF SO THEN
!  TAKE APPROPRIATE ACTION
!
      DO II=1,L
        DO JI=II+1,L
          IF(ABS(A(II,JI)).GT.TOL2)IOFF=IOFF+1
          IF(ABS(A(JI,II)).GT.TOL2)IOFF=IOFF+1
        END DO
      END DO
      IF(IOFF.EQ.0)THEN
          CALL ZERO(B,L*L)
          DO 40 I=1,L
          B(I,I)=ONE
40        CONTINUE
      ELSE
50    P=ZER
      DO 70 I=2,N
      IM1=I-1
      DO 60 J=1,IM1
      Q=A(I,J)
      IF(P.GE. ABS(Q)) GO TO 60
      P= ABS(Q)
      II=I
      JJ=J
60    CONTINUE
70    CONTINUE
      IF(P.EQ.0.) GO TO 140
      P=A(II,II)
      Q=A(II,JJ)
      R=A(JJ,JJ)
      DIFF=0.5D0*(P-R)
      IF( ABS(DIFF).LT. ABS(Q)) GO TO 80
      IF( ABS(Q/DIFF).GT.TOL) GO TO 80
      A(II,JJ)=ZER
      A(JJ,II)=ZER
      GO TO 50
80    S=SQRT(0.250D0*(P-R)**2+Q**2)
      SUM=0.5D0*(P+R)
      D=R*P-Q**2
      IF(SUM.GT.ZER) GO TO 90
      ALN=SUM-S
      ALP=D/ALN
      GO TO 100
90    ALP=SUM+S
      ALN=D/ALP
100   IF(DIFF.GT.ZER) GO TO 110
      T=Q/(DIFF-S)
      A(II,II)=ALN
      A(JJ,JJ)=ALP
      GO TO 120
110   T=Q/(DIFF+S)
      A(II,II)=ALP
      A(JJ,JJ)=ALN
120   C=1.0D0/SQRT(1.0D0+T**2)
      S=T*C
      A(II,JJ)=ZER
      A(JJ,II)=ZER
      DO 130 I=1,N
      P=B(I,II)
      Q=B(I,JJ)
      B(I,II)=C*P+S*Q
      B(I,JJ)=C*Q-S*P
      IF(I.EQ.II.OR.I.EQ.JJ) GO TO 130
      P=A(I,II)
      Q=A(I,JJ)
      R=C*P+S*Q
       A(I,II)=R
      A(II,I)=R
      R=Q*C-P*S
      A(I,JJ)=R
      A(JJ,I)=R
130   CONTINUE
      GO TO 50
      ENDIF
140   IF(N1.EQ.1) RETURN
      MM=N-1
      DO I=1,MM
        II=I+1
        DO J=II,N
          IF(A(I,I)-A(J,J) .LT. 0) THEN
            CYCLE
          ELSE 
            GOTO 150
          END IF
150       W1=A(I,I)
          A(I,I)=A(J,J)
          A(J,J)=W1
          DO MU=1,N
            W2=B(MU,I)
            B(MU,I)=B(MU,J)
            B(MU,J)=W2
          END DO
        END DO
      END DO
      RETURN
      END



      SUBROUTINE ZERO(A,NA)
      IMPLICIT NONE
      INTEGER NA,I
      DOUBLE PRECISION A(NA)

      DO I=1,NA
        A(I)=0.D0
      END DO
      CONTINUE
      RETURN
      END

