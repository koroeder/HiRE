!  Subroutine CENTRE moves the centre of mass to the origin.

      SUBROUTINE CENTRE2(X)
      USE PREC
      USE COMMONS, ONLY: NATOMS
      USE GENRIGID
      IMPLICIT NONE
      REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
      REAL(KIND = REAL64) :: XMASS, YMASS, ZMASS
      INTEGER I, J1
         
      ! if generalised rigid body is used, be careful when averaging
      ! no need to shift the rotational degrees of freedom!
      IF ( .NOT. ATOMRIGIDCOORDT ) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NRIGIDBODY
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               XMASS=XMASS+X(6*NRIGIDBODY + 3*J1-2)
               YMASS=YMASS+X(6*NRIGIDBODY + 3*J1-1)
               ZMASS=ZMASS+X(6*NRIGIDBODY + 3*J1  )
            ENDDO
         ENDIF
         XMASS=XMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         YMASS=YMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         ZMASS=ZMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         DO I=1,NRIGIDBODY
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               X(6*NRIGIDBODY + 3*J1-2) = X(6*NRIGIDBODY + 3*J1-2) - XMASS
               X(6*NRIGIDBODY + 3*J1-1) = X(6*NRIGIDBODY + 3*J1-1) - YMASS
               X(6*NRIGIDBODY + 3*J1) = X(6*NRIGIDBODY+3*J1) - ZMASS
            ENDDO
         ENDIF

      ELSE
      !proceeds as usual if everything is in atom coords
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
         DO I=1,NATOMS
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
      ENDIF

      RETURN
      END SUBROUTINE CENTRE2

      SUBROUTINE SETCENTRE(X)
      USE PREC
      USE COMMONS, ONLY: NATOMS, DEBUG, MYUNIT, CENTX, CENTY, CENTZ
      USE GENRIGID

      IMPLICIT NONE
      REAL(KIND = REAL64), INTENT(INOUT) :: X(3*NATOMS)
      REAL(KIND = REAL64) XMASS, YMASS, ZMASS
      INTEGER I, J1

      IF ( .NOT. ATOMRIGIDCOORDT ) THEN
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NRIGIDBODY
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               XMASS=XMASS+X(6*NRIGIDBODY + 3*J1-2)
               YMASS=YMASS+X(6*NRIGIDBODY + 3*J1-1)
               ZMASS=ZMASS+X(6*NRIGIDBODY + 3*J1  )
            ENDDO
         ENDIF
         XMASS=XMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         YMASS=YMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         ZMASS=ZMASS/ ( NRIGIDBODY + (DEGFREEDOMS - 6*NRIGIDBODY)/3)
         DO I=1,NRIGIDBODY
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS
         ENDDO
         IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
            DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
               X(6*NRIGIDBODY + 3*J1-2) = X(6*NRIGIDBODY + 3*J1-2) - XMASS
               X(6*NRIGIDBODY + 3*J1-1) = X(6*NRIGIDBODY + 3*J1-1) - YMASS
               X(6*NRIGIDBODY+3*J1) = X(6*NRIGIDBODY+3*J1) - ZMASS
            ENDDO
         ENDIF

      ELSE
         !XMASS, YMASS and ZMASS are the components of the COM position vector
         XMASS=0.0D0
         YMASS=0.0D0
         ZMASS=0.0D0
         DO I=1,NATOMS
            XMASS=XMASS+X(3*(I-1)+1)
            YMASS=YMASS+X(3*(I-1)+2)
            ZMASS=ZMASS+X(3*(I-1)+3)
         ENDDO
         XMASS=XMASS/NATOMS
         YMASS=YMASS/NATOMS
         ZMASS=ZMASS/NATOMS
         DO I=1,NATOMS
            X(3*(I-1)+1)=X(3*(I-1)+1)-XMASS+CENTX
            X(3*(I-1)+2)=X(3*(I-1)+2)-YMASS+CENTY
            X(3*(I-1)+3)=X(3*(I-1)+3)-ZMASS+CENTZ
         ENDDO
         IF (DEBUG) WRITE(MYUNIT,'(A,3F12.4)') 'centre of mass moved to ',CENTX,CENTY,CENTZ
      ENDIF

      RETURN
      END SUBROUTINE SETCENTRE

