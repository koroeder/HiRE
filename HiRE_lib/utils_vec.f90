module vec_utils
use prec_hire
contains

  ! Euclidean norm of a vector
  pure function euc_norm(v)
    real(kind = real64) :: euc_norm
    real(kind = real64), intent(in) :: v(3)

    euc_norm = dsqrt(dot_product(v,v))
  end function euc_norm
  
  ! Cross product
  pure function crossproduct(u, v)
    real(kind = real64) :: crossproduct(3)
    real(kind = real64), intent (in) :: u(3), v(3)
  
    crossproduct(1) = u(2)*v(3) - u(3)*v(2)
    crossproduct(2) = u(3)*v(1) - u(1)*v(3)
    crossproduct(3) = u(1)*v(2) - u(2)*v(1) 
  end function crossproduct

  ! Normalised cross product
  subroutine normed_cp(u, v, uxv, norm)
    real(kind = real64), intent(out) :: uxv(3)
    real(kind = real64), intent(in)  :: u(3), v(3)  
    real(kind = real64), intent(out) :: norm
    
    uxv(1:3) = crossproduct(u,v)
    norm = euc_norm(uxv)
    uxv(1:3) = uxv(1:3)/norm
  end subroutine normed_cp

  ! Normalised crossproduct, also returning unnormalised cp
  subroutine normed_cp2(u, v, uxv, uxv0, norm)
    real(kind = real64), intent(out) :: uxv(3), uxv0(3)
    real(kind = real64), intent(in)  :: u(3), v(3)  
    real(kind = real64), intent(out) :: norm
    
    uxv(1:3) = crossproduct(u,v)
    norm = euc_norm(uxv)
    uxv0(1:3) = uxv(1:3)/norm
  end subroutine normed_cp2


  ! subroutine to get norm and normed vector
  subroutine normed_vec(vec, vec0, norm)
    real(kind = real64), intent(in) :: vec(3)
    real(kind = real64), intent(out) :: vec0(3)  
    real(kind = real64), intent(out) :: norm
    
    norm = euc_norm(vec)
    vec0(1:3) = vec(1:3)/norm   
  end subroutine normed_vec

  ! difference between two vectors
  pure function vec_diff(v1, v2) result(v3)
     implicit none
     real(kind = real64), intent(in) :: v1(3), v2(3)
     real(kind = real64) v3(3)
     v3(1) = v2(1) - v1(1)
     v3(2) = v2(2) - v1(2)
     v3(3) = v2(3) - v1(3)
  end function vec_diff
  
end module vec_utils
