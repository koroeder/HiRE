subroutine calcforce_RNA_protein(scale, X, F, E_cross)
  use definitions, only: vecsize, Natoms, N_prot
  use prec_hire
  use geometric_corrections
  implicit none
  real(kind = real64), intent(in) :: scale
  real(kind = real64), dimension(vecsize), intent(in) :: X
  real(kind = real64), dimension(vecsize), intent(out) :: F
  real(kind = real64), intent(out) :: E_cross

  integer :: i, j
  real(kind = real64) :: r(3), d2, dm6, E_c, F_c

  do i=1, N_prot
    do j=N_prot+1, Natoms
      r = X(i*3-2:i*3) - X(j*3-2:j*3)
      d2 = dot_product(r,r)
      dm6 = 1/d2**3
      ! simple lj for now
      E_c = dm6*(dm6-1)
      E_cross = E_cross + E_c
      F_c = 6*(dm6*(2*dm6-1))/d2
      F(i*3-2:i*3) = F(i*3-2:i*3) + F_c*r 
      F(j*3-2:j*3) = F(j*3-2:j*3) - F_c*r 
    enddo
  enddo

end subroutine calcforce_RNA_protein

