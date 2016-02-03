module derivatives
  use numbers
contains
  subroutine constant_accel(x,y,dydx)
    implicit none
    real(dp), intent(in)  :: x, y(:)
    real(dp), intent(out) :: dydx(:)
    dydx = -9.81
  end subroutine constant_accel
end module derivatives
