module dyn_sys
  use nrtype
  implicit none
  real(dp) :: a, b, c, u, d, alpha, beta

contains

  subroutine Chen(x,y,dydx)
    ! https://en.wikipedia.org/wiki/Multiscroll_attractor
    ! x := t
    ! y(1:3) := x, y, z
    ! dydx(1:3) := dx/dt, dy/dt, dz/dt
    implicit none
    real(dp), intent(in)  :: x, y(:)
    real(dp), intent(out) :: dydx(:)

    dydx(1) = a * (y(2) - y(1))
    dydx(2) = (c-a) * y(1) - y(1) * y(3) + c * y(2)
    dydx(3) = y(1) * y(2) - b * y(3)
  end subroutine Chen

  subroutine Lu_Chen(x,y,dydx)
    ! https://en.wikipedia.org/wiki/Multiscroll_attractor
    ! x := t
    ! y(1:3) := x, y, z
    ! dydx(1:3) := dx/dt, dy/dt, dz/dt
    implicit none
    real(dp), intent(in)  :: x, y(:)
    real(dp), intent(out) :: dydx(:)

    dydx(1) = a * (y(2) - y(1))
    dydx(2) = y(1) - y(1) * y(3) + c * y(2) + u
    dydx(3) = y(1) * y(2) - b * y(3)
  end subroutine Lu_Chen

  subroutine mc_Chua(x,y,dydx)
    implicit none
    real(dp), intent(in)  :: x, y(:)
    real(dp), intent(out) :: dydx(:)

    dydx(1) = alpha * (y(2) - b * sin(pi*y(1))/(2.*a) + d)
    dydx(2) = y(1) - y(2) + y(3)
    dydx(3) = -beta * y(2)
  end subroutine mc_Chua

end module dyn_sys
