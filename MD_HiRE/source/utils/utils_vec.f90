!> @file
!> File contains UTILS_VEC module with basic vector routines

!> Module containing routines to deal with vectors
!> @brief

MODULE UTILS_VEC
   USE NUMKIND
   IMPLICIT NONE

   CONTAINS
      !> Function to get vector length
      !> @brief
      PURE FUNCTION VEC_LEN(V)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: V(3)
         REAL(KIND=REAL64) :: VEC_LEN

         VEC_LEN = SQRT(DOT_PRODUCT(V,V))
      END FUNCTION VEC_LEN

      !> Function to get normalised vector
      PURE FUNCTION VEC_NORMED(V) RESULT(U)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: V(3)
         REAL(KIND=REAL64) :: U(3), NORM

         NORM = VEC_LEN(V)
         U(1:3) = V(1:3)/NORM
      END FUNCTION VEC_NORMED
         

      !> Function to get difference between vectors
      !> @brief
      !>
      !> Returns difference between V1 and V2 as V2-V1
      !>
      !> @param[in] V1 - vector 1
      !> @param[in] V2 - vector 2
      !>
      !> @return Returns difference between vectors
      PURE FUNCTION VEC_DIFF(V1,V2) RESULT(V3)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: V1(3), V2(3)
         REAL(KIND=REAL64) :: V3(3)

         V3(1) = V2(1) - V1(1)
         V3(2) = V2(2) - V1(2)
         V3(3) = V2(3) - V1(3)
      END FUNCTION VEC_DIFF

      !> Function to get crossproduct of vectors
      !> @brief
      !>
      !> Returns crossproduct between U and V as U x V
      !>
      !> @param[in] V - vector 1
      !> @param[in] U - vector 2
      !>
      !> @return Returns the crossproduct of vectors     
      PURE FUNCTION CROSSP(U, V)
         IMPLICIT NONE
         REAL(KIND=REAL64), INTENT(IN) :: U(3), V(3)
         REAL(KIND=REAL64) :: CROSSP(3)

         CROSSP(1) = U(2)*V(3) - U(3)*V(2)
         CROSSP(2) = U(3)*V(1) - U(1)*V(3)
         CROSSP(3) = U(1)*V(2) - U(2)*V(1)
      END FUNCTION

      !> convert angle axis to quaternion
      function rot_aa2q(p) result(q)
         implicit none

         real(kind = real64), intent(in) :: p(3)
         real(kind = real64) :: q(4)
         real(kind = real64) :: thetah
         real(kind = real64), parameter :: epsilon = 1.0d-6

         thetah = 0.5d0 * vec_len(p)
         q(1) = cos(thetah)

         ! do linear expansion for small epsilon
         if (thetah < epsilon) then
            q(2:4) = 0.5d0 * p
         else
            q(2:4) = 0.5d0 * sin(thetah) * p / thetah
         endif
         ! make sure to have normal form
         if(q(1) < 0d0) q = -q
      end function

      !> convert quarternion into rot matrix
      function rot_q2mx(qin) result(m)
         implicit none
         real(kind = real64), intent(in) :: qin(4)
         real(kind = real64) :: m(3,3)
         real(kind = real64) :: q(4)
         real(kind = real64) :: sq(4), tmp1, tmp2
         integer i

         q = qin / sqrt(dot_product(qin,qin))

         do i=1,4
            sq(i) = q(i)*q(i)
         enddo

         m(1,1) = ( sq(2) - sq(3) - sq(4) + sq(1))
         m(2,2) = (-sq(2) + sq(3) - sq(4) + sq(1))
         m(3,3) = (-sq(2) - sq(3) + sq(4) + sq(1))

         tmp1 = q(2)*q(3)
         tmp2 = q(1)*q(4)
         m(2,1) = 2.0d0 * (tmp1 + tmp2)
         m(1,2) = 2.0d0 * (tmp1 - tmp2)

         tmp1 = q(2)*q(4)
         tmp2 = q(3)*q(1)
         m(3,1) = 2.0d0 * (tmp1 - tmp2)
         m(1,3) = 2.0d0 * (tmp1 + tmp2)
         tmp1 = q(3)*q(4)
         tmp2 = q(1)*q(2)
         m(3,2) = 2.0d0 * (tmp1 + tmp2)
         m(2,3) = 2.0d0 * (tmp1 - tmp2)
      end function
      
      !> convert aa to rot matrix
      function rot_aa2mx(p) result(m)
         implicit none
         real(kind = real64), intent(in) :: p(3)
         real(kind = real64) m(3,3)

         m=rot_q2mx(rot_aa2q(p))
      end function

END MODULE UTILS_VEC