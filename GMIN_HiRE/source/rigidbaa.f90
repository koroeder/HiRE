! Finds the rotation matrix and its derivatives with respect to 
! components of the rotation vector.
!
! See 'Simulations of rigid bodies in an angle axis framework', Dwaipayan
! Chakrabarti and David Wales, Phys. Chem. Chem. Phys., 2009, 11, 
! 1970-1976
! Modified for general angle-axis 30/01/12
! 

      SUBROUTINE ROTMAT(P, RM)

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, CT, ST, I3(3,3), E(3,3), RM(3,3)

      I3(:,:) = 0.D0
      DO K = 1, 3
         I3(K,K) = 1.D0
      ENDDO

      THETA2 = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN
         RM(:,:) = I3(:,:)
      ELSE
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         RM      = I3(:,:) + (1.D0-CT)*MATMUL(E(:,:),E(:,:)) + ST*E(:,:)
      ENDIF

      END SUBROUTINE ROTMAT

!     --------------------------------------------------------------------------

!     RMDVDT = rotation matrix derivative
!     P(3) = rotation vector 
!     RM(3,3) = rotation matrix
!     DRMk(3,3) = derivative of the rotation matrix with respect to the 
!                 kth component of the rotation vector
!     GTEST = true if derivatives are to be found  

      SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)

      IMPLICIT NONE

!     PN(3) = the unit vector parallel to P
!     THETA = the modulus of the rotation vector P, equivalent to the 
!             angle of rotation
!     THETA2 = THETA squared
!     THETA3 = THETA**-3
!     CT = cos(THETA)
!     ST = sin(THETA)
!     I3(3,3) = 3x3 identity matrix
!     E(3,3) = the skew-symmetric matrix obtained from a unit vector 
!              parallel to P (equation (2) in the paper)
!     ESQ(3,3) = the square of E
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)

!     DEk = derivate of E with respect to the kth component of P
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      LOGICAL, intent(in) :: GTEST

!     Set the values of the idenity matrix I3
      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!     Calculate the value of THETA2 as the square modulus of P
      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 < 1.0D-12) THEN
!        Execute if the angle of rotation is zero
!        In this case the rotation matrix is the identity matrix
         RM(:,:) = I3(:,:)

         ! vr274> first order corrections to rotation matrix
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

!        If derivatives do not need to found, we're finished
         IF (.NOT. GTEST) RETURN

!        This is the special case described in the paper, where DRMk is
!        equal to E(k), which is the skew-symmetric matrix obtained from
!        P with P(k) equal to 1 and other components equal to zero
!         PN        = (/1.D0, 0.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(2,3)    = -PN(1)
!         E(3,2)    = -E(2,3)
!         DRM1(:,:) = E(:,:)

!         PN        = (/0.D0, 1.D0, 0.D0/)
!         E(:,:)    = 0.D0
!         E(1,3)    =  PN(2)
!         E(3,1)    = -E(1,3)
!         DRM2(:,:) = E(:,:)

!         PN        = (/0.D0, 0.D0, 1.D0/)
!         E(:,:)    = 0.D0
!         E(1,2)    = -PN(3)
!         E(2,1)    = -E(1,2)
!         DRM3(:,:) = E(:,:)

! now up to the linear order in theta
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    = -2.0D0*P(1)
         E(2,3)    = -2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    = -2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    = -2.0D0
         E(3,2)    = P(3)
         E(3,3)    = -2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(3)
         E(1,2)    = -2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    = -2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

      ELSE
!       Execute for the general case, where THETA dos not equal zero
!       Find values of THETA, CT, ST and THETA3
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)

!        Set THETA to 1/THETA purely for convenience
         THETA   = 1.D0/THETA

!        Normalise P and construct the skew-symmetric matrix E
!        ESQ is calculated as the square of E
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         ESQ     = MATMUL(E,E)

!        RM is calculated from Rodrigues' rotation formula (equation (1)
!        in the paper)
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

!        If derivatives do not need to found, we're finished
         IF (.NOT. GTEST) RETURN

!        Set up DEk using the form given in equation (4) in the paper
         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

!        Use equation (3) in the paper to find DRMk
         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

      ENDIF
      END SUBROUTINE RMDRVT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RMDFAS(P, RM, DRM1, DRM2, DRM3, D2RM1, D2RM2, D2RM3, D2RI12, D2RI23, D2RI31, GTEST, STEST)

      IMPLICIT NONE

      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, E(3,3), ESQ(3,3), I3(3,3)
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      DOUBLE PRECISION :: D2E1(3,3), D2E2(3,3), D2E3(3,3), D2E12(3,3), D2E23(3,3), D2E31(3,3)
      DOUBLE PRECISION :: D2RM1(3,3), D2RM2(3,3), D2RM3(3,3)
      DOUBLE PRECISION :: D2RI12(3,3), D2RI23(3,3), D2RI31(3,3)
      DOUBLE PRECISION :: FCTR, FCTRSQ1, FCTRSQ2, FCTRSQ3
      LOGICAL          :: GTEST, STEST

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!      RM(3,3) = 0.D0; DRM1(3,3) = 0.D0; DRM2(3,3) = 0.D0; DRM3(3,3) = 0.D0
!      D2RM1(3,3) = 0.D0; D2RM2(3,3) = 0.D0; D2RM3(3,3) = 0.D0
!      D2RI12(3,3) = 0.D0; D2RI23(3,3) = 0.D0; D2RI31(3,3) = 0.D0

      THETA2  = DOT_PRODUCT(P,P)

      IF (ABS(THETA2) < 1.D-12) THEN

         RM(:,:)     = I3(:,:)
! first order corrections to rotation matrix
! 11/1/12
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

! now up to the linear order in theta
! 11/1/12
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    = -2.0D0*P(1)
         E(2,3)    = -2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    = -2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)
         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    = -2.0D0
         E(3,2)    = P(3)
         E(3,3)    = -2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(3)
         E(1,2)    = -2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    = -2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

         IF (.NOT. STEST) RETURN

! Fix bug and now to linear order in theta
! 11/1/12
         D2RM1(:,:) = 0.0D0
         D2RM1(1,1) = 0.0D0
         D2RM1(1,2) =  P(3)
         D2RM1(1,3) = -P(2)
         D2RM1(2,1) = -P(3)
         D2RM1(2,2) = -3.0D0
         D2RM1(2,3) =  3.0D0*P(1)
         D2RM1(3,1) =  P(2)
         D2RM1(3,2) = -3.0D0*P(1)
         D2RM1(3,3) = -3.0D0
         D2RM1(:,:) = D2RM1(:,:)/3.0D0

         D2RM2(:,:) = 0.0D0
         D2RM2(1,1) = -3.0D0
         D2RM2(1,2) =  P(3)
         D2RM2(1,3) = -3.0D0*P(2)
         D2RM2(2,1) = -P(3)
         D2RM2(2,2) =  0.0D0
         D2RM2(2,3) =  P(1)
         D2RM2(3,1) =  3.0D0*P(2)
         D2RM2(3,2) = -P(1)
         D2RM2(3,3) = -3.0D0
         D2RM2(:,:) = D2RM2(:,:)/3.0D0

        D2RM3(:,:) = 0.0D0
         D2RM3(1,1) = -3.0D0
         D2RM3(1,2) =  3.0D0*P(3)
         D2RM3(1,3) = -P(2)
         D2RM3(2,1) = -3.0D0*P(3)
         D2RM3(2,2) = -3.0D0
         D2RM3(2,3) =  P(1)
         D2RM3(3,1) =  P(2)
         D2RM3(3,2) = -P(1)
         D2RM3(3,3) =  0.0D0
         D2RM3(:,:) = D2RM3(:,:)/3.0D0

         D2RI12(:,:) = 0.0D0
         D2RI12(1,1) = 0.0D0
         D2RI12(1,2) = 0.5D0
         D2RI12(1,3) = -P(1)/3.0D0
         D2RI12(2,1) = 0.5D0
         D2RI12(2,2) = 0.0D0
         D2RI12(2,3) = P(2)/3.0D0
         D2RI12(3,1) = P(1)/3.0D0
         D2RI12(3,2) = -P(2)/3.0D0
         D2RI12(3,3) = 0.0D0

         D2RI23(:,:) = 0.0D0
         D2RI23(1,1) = 0.0D0
         D2RI23(1,2) = P(2)/3.0D0
         D2RI23(1,3) = -P(3)/3.0D0
         D2RI23(2,1) = -P(2)/3.0D0
        D2RI23(2,2) = 0.0D0
         D2RI23(2,3) = 0.5D0
         D2RI23(3,1) = P(3)/3.0D0
         D2RI23(3,2) = 0.5D0
         D2RI23(3,3) = 0.0D0

         D2RI31(:,:) = 0.0D0
         D2RI31(1,1) = 0.0D0
         D2RI31(1,2) = P(1)/3.0D0
         D2RI31(1,3) = 0.5D0
         D2RI31(2,1) = -P(1)/3.0D0
         D2RI31(2,2) = 0.0D0
         D2RI31(2,3) = P(3)/3.0D0
         D2RI31(3,1) = 0.5D0
         D2RI31(3,2) = -P(3)/3.0D0
         D2RI31(3,3) = 0.0D0

      ELSE

         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)
        THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)

         ESQ     = MATMUL(E(:,:),E(:,:))
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
        DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

         IF (.NOT. STEST) RETURN

         FCTRSQ1 = PN(1)*PN(1)
         FCTRSQ2 = PN(2)*PN(2)
         FCTRSQ3 = PN(3)*PN(3)
        D2E1(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(1)*PN(1))*THETA3
         D2E1(1,2) =  P(3)*FCTR
         D2E1(1,3) = -P(2)*FCTR
         D2E1(2,3) =  P(1)*FCTR + 2.D0*P(1)*THETA3
         D2E1(2,1) = -D2E1(1,2)
         D2E1(3,1) = -D2E1(1,3)
         D2E1(3,2) = -D2E1(2,3)

         D2E2(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(2)*PN(2))*THETA3
         D2E2(1,2) =  P(3)*FCTR
         D2E2(1,3) = -P(2)*FCTR - 2.D0*P(2)*THETA3
         D2E2(2,3) =  P(1)*FCTR
         D2E2(2,1) = -D2E2(1,2)
         D2E2(3,1) = -D2E2(1,3)
         D2E2(3,2) = -D2E2(2,3)

         D2E3(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(3)*PN(3))*THETA3
         D2E3(1,2) =  P(3)*FCTR + 2.D0*P(3)*THETA3
         D2E3(1,3) = -P(2)*FCTR
         D2E3(2,3) =  P(1)*FCTR
         D2E3(2,1) = -D2E3(1,2)
         D2E3(3,1) = -D2E3(1,3)
         D2E3(3,2) = -D2E3(2,3)
         D2E12(:,:) =  0.D0
         D2E12(1,2) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E12(1,3) = -PN(1)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E12(2,3) =  PN(2)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E12(2,1) = -D2E12(1,2)
         D2E12(3,1) = -D2E12(1,3)
         D2E12(3,2) = -D2E12(2,3)

         D2E23(:,:) =  0.D0
         D2E23(1,2) =  PN(2)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E23(1,3) = -PN(3)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E23(2,3) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E23(2,1) = -D2E23(1,2)
         D2E23(3,1) = -D2E23(1,3)
         D2E23(3,2) = -D2E23(2,3)

         D2E31(:,:) =  0.D0
         D2E31(1,2) =  PN(1)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E31(1,3) =  3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E31(2,3) =  PN(3)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E31(2,1) = -D2E31(1,2)
         D2E31(3,1) = -D2E31(1,3)
         D2E31(3,2) = -D2E31(2,3)

         D2RM1(:,:)  = ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (CT*FCTRSQ1 - ST*FCTRSQ1*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE1(:,:),DE1(:,:)) + MATMUL(D2E1(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E1(:,:))) + (- ST*FCTRSQ1 - CT*FCTRSQ1*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(1)*DE1(:,:) + ST*D2E1(:,:)

         D2RM2(:,:)  = ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (CT*FCTRSQ2 - ST*FCTRSQ2*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE2(:,:),DE2(:,:)) + MATMUL(D2E2(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E2(:,:))) + (- ST*FCTRSQ2 - CT*FCTRSQ2*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(2)*DE2(:,:) + ST*D2E2(:,:)

         D2RM3(:,:)  = ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (CT*FCTRSQ3 - ST*FCTRSQ3*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (1.D0-CT)*(2.D0*MATMUL(DE3(:,:),DE3(:,:)) + MATMUL(D2E3(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E3(:,:))) + (- ST*FCTRSQ3 - CT*FCTRSQ3*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(3)*DE3(:,:) + ST*D2E3(:,:)

         D2RI12(:,:) = ST*PN(1)*(MATMUL(E(:,:),DE2(:,:)) + MATMUL(DE2(:,:),E(:,:))) &
                     + (CT*PN(1)*PN(2) - ST*PN(1)*PN(2)*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E12(:,:),E(:,:)) + MATMUL(DE1(:,:),DE2(:,:)) &
                     + MATMUL(DE2(:,:),DE1(:,:)) + MATMUL(E(:,:),D2E12(:,:))) &
                     - (ST*PN(1)*PN(2) + CT*PN(1)*PN(2)*THETA)*E(:,:) + CT*PN(1)*DE2(:,:) + CT*PN(2)*DE1(:,:) &
                     + ST*D2E12(:,:)

         D2RI23(:,:) = ST*PN(2)*(MATMUL(E(:,:),DE3(:,:)) + MATMUL(DE3(:,:),E(:,:))) &
                     + (CT*PN(2)*PN(3) - ST*PN(2)*PN(3)*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E23(:,:),E(:,:)) + MATMUL(DE2(:,:),DE3(:,:)) &
                     + MATMUL(DE3(:,:),DE2(:,:)) + MATMUL(E(:,:),D2E23(:,:))) &
                     - (ST*PN(2)*PN(3) + CT*PN(2)*PN(3)*THETA)*E(:,:) + CT*PN(2)*DE3(:,:) + CT*PN(3)*DE2(:,:) &
                     + ST*D2E23(:,:)

         D2RI31(:,:) = ST*PN(3)*(MATMUL(E(:,:),DE1(:,:)) + MATMUL(DE1(:,:),E(:,:))) &
                     + (CT*PN(3)*PN(1) - ST*PN(3)*PN(1)*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E31(:,:),E(:,:)) + MATMUL(DE3(:,:),DE1(:,:)) &
                     + MATMUL(DE1(:,:),DE3(:,:)) + MATMUL(E(:,:),D2E31(:,:))) &
                     - (ST*PN(3)*PN(1) + CT*PN(3)*PN(1)*THETA)*E(:,:) + CT*PN(3)*DE1(:,:) + CT*PN(1)*DE3(:,:) &
                     + ST*D2E31(:,:)

      ENDIF

      END SUBROUTINE RMDFAS

!     ----------------------------------------------------------------------------------------------


      SUBROUTINE QROTMAT(Q,RM)

!     provides the rotation matrix RM corresponding to quaternion Q
!     right-handed rotation in the right-handed coordinate system

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), RM(3,3)

      RM(1,1) = Q(1)**2 + Q(2)**2 - Q(3)**2 - Q(4)**2
      RM(1,2) = 2.D0*(Q(2)*Q(3) - Q(1)*Q(4))
      RM(1,3) = 2.D0*(Q(2)*Q(4) + Q(1)*Q(3))
      RM(2,1) = 2.D0*(Q(2)*Q(3) + Q(1)*Q(4))
      RM(2,2) = Q(1)**2 + Q(3)**2 - Q(2)**2 - Q(4)**2
      RM(2,3) = 2.D0*(Q(3)*Q(4) - Q(1)*Q(2))
      RM(3,1) = 2.D0*(Q(2)*Q(4) - Q(1)*Q(3))
      RM(3,2) = 2.D0*(Q(3)*Q(4) + Q(1)*Q(2))
      RM(3,3) = Q(1)**2 + Q(4)**2 - Q(2)**2 - Q(3)**2

      END SUBROUTINE





