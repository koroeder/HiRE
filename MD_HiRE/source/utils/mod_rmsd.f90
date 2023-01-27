MODULE MOD_RMSD
   USE NUMKIND
   IMPLICIT NONE
   REAL(KIND = REAL64), ALLOCATABLE :: REFX(:)
   REAL(KIND = REAL64) :: REFCX(3)
   CONTAINS
      
      SUBROUTINE SET_REF(NATOMS,X)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)

         IF (.NOT.ALLOCATED(REFX)) ALLOCATE(REFX(3*NATOMS))
         ! set reference to new coordinates
         REFX(1:3*NATOMS) = X(1:3*NATOMS)
         ! find origin and centre the coordinates
         CALL FIND_ORIGIN(NATOMS,REFX,REFCX)
         CALL CENTRE_COORDS(NATOMS,REFX,REFCX)
      END SUBROUTINE SET_REF

      SUBROUTINE GET_RMSD(NATOMS, X, DIST, RMSD, ALIGNT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)      
         REAL(KIND = REAL64), INTENT(OUT) :: DIST, RMSD
         LOGICAL, INTENT(IN) :: ALIGNT
         REAL(KIND = REAL64) :: RMAT(3,3) 

         CALL FIND_ALIGNMENT(NATOMS, X, DIST, RMAT)
         RMSD = DIST/SQRT(DBLE(NATOMS))

         IF (ALIGNT) THEN
            CALL ALIGNTOREF(NATOMS, X, RMAT)
         END IF

      END SUBROUTINE GET_RMSD

      SUBROUTINE FIND_ALIGNMENT(NATOMS, X, DIST, RMAT)
         USE MD_COMMONS, ONLY: MYUNIT
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS) 
         REAL(KIND = REAL64), INTENT(OUT) :: DIST
         REAL(KIND = REAL64), INTENT(OUT) :: RMAT(3,3)   
         INTEGER, PARAMETER :: LWORK=12       
         REAL(KIND = REAL64) :: CX(3)
         REAL(KIND = REAL64) :: QMAT(4,4) , XM, YM, ZM, XP, YP, ZP, DIAG(4), TEMPA(LWORK), MINV
         REAL(KIND = REAL64) :: Q1, Q2, Q3, Q4
         INTEGER :: I, J, JMIN, INFO
         ! move coordinates to the centre
         CALL FIND_ORIGIN(NATOMS,X,CX)
         CALL CENTRE_COORDS(NATOMS,X,CX)
         ! Analystic method based on quaternions and general angle-axis
         ! See: Kearsley, Acta Cryst. A, 45, 208-210, 1989
         !      Griffiths, Niblett and Wales, JCTC, 13, 4914-1931, 2017
         QMAT(1:4,1:4)=0.0D0
         DO I=1,NATOMS
            XM=REFX(3*(I-1)+1)-X(3*(I-1)+1)
            YM=REFX(3*(I-1)+2)-X(3*(I-1)+2)
            ZM=REFX(3*(I-1)+3)-X(3*(I-1)+3)
            XP=REFX(3*(I-1)+1)+X(3*(I-1)+1)
            YP=REFX(3*(I-1)+2)+X(3*(I-1)+2)
            ZP=REFX(3*(I-1)+3)+X(3*(I-1)+3)
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
         QMAT(2,1)=QMAT(1,2)
         QMAT(3,1)=QMAT(1,3)
         QMAT(3,2)=QMAT(2,3)
         QMAT(4,1)=QMAT(1,4)
         QMAT(4,2)=QMAT(2,4)
         QMAT(4,3)=QMAT(3,4)
         !Eigendecomposition fo the quarternion
         CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,LWORK,INFO)
         MINV=1.0D100
         DO J=1,4
            IF (DIAG(J).LT.MINV) THEN
               JMIN=J
               MINV=DIAG(J)
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
         ! This is the Euclidean distance!
         DIST=SQRT(MINV)
         ! Get the rotational matrix
         Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
         RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
         RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
         RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
         RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
         RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
         RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
         RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
         RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
         RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2             
      END SUBROUTINE FIND_ALIGNMENT

      SUBROUTINE ALIGNTOREF(NATOMS,X,RMAT)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
         REAL(KIND = REAL64), INTENT(IN) :: RMAT(3,3)
         INTEGER :: I, J, K
         REAL(KIND = REAL64) :: R1, R0(3)

         ! translate the coordinates to reference centre, 
         ! the structure itself was shifted to the origin in the FIND_ALIGNMENT routine
         CALL CENTRE_COORDS(NATOMS,X,-REFCX)

         DO I=1,NATOMS
            DO J=1,3
               R0(J) = X(3*(I-1)+1) - REFCX(J)
            END DO
            DO J=1,3
               R1=0.0D0
               DO K=1,3
                  R1 = R1 + RMAT(J,K)*R0(K)
               ENDDO
               X(3*(I-1)+J) = R1 + REFCX(J)
            ENDDO
         ENDDO     
      END SUBROUTINE ALIGNTOREF
        
      SUBROUTINE FIND_ORIGIN(NATOMS,X,CX)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(IN) :: X(3*NATOMS)
         REAL(KIND = REAL64), INTENT(OUT) :: CX(3)
         INTEGER :: I, J

         CX(1:3) = 0.0D0

         DO I=1,NATOMS
            DO J = 1,3
               CX(J) = CX(J) + X(3*(I-1) + J)
            END DO
         END DO
         CX(1:3) = CX(1:3)/DBLE(NATOMS)
      END SUBROUTINE FIND_ORIGIN

      SUBROUTINE CENTRE_COORDS(NATOMS,X,CX)
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
         REAL(KIND = REAL64), INTENT(IN) :: CX(3)
         INTEGER :: I, J, IDX
         
         DO I=1,NATOMS
            DO J=1,3
               IDX = 3*(I-1)+J
               X(IDX) = X(IDX) - CX(J)
            END DO
         END DO
      END SUBROUTINE CENTRE_COORDS

END MODULE MOD_RMSD