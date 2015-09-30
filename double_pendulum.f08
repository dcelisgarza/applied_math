module d_pend
  double precision :: l1, l2, m1, m2, ms
end module d_pend

program main
  use d_pend
  implicit none
  call d_pend_stat
  call d_pend_adapt
end program main

subroutine d_pend_stat
  use d_pend
  use ode_int
  implicit none
  integer, parameter :: n = 4
  double precision, parameter :: rad = atan(1.)/45.
  double precision t, dt, tmax, dtheta
  double precision yi(n)
  double precision, dimension(size(yi)) :: yf, dydt, er
  double precision start, finish
  interface
    subroutine d_pend_deriv(t,y,dydt)
      implicit none
      double precision, intent(in) :: t, y(:)
      double precision, intent(out) :: dydt(:)
    end subroutine d_pend_deriv
  end interface

  open(unit=1, file = 'd_pend_stat.plt')
  ! Initial conditions.
  t = 0.
  tmax = 50.
  dt = 5.e-3
  l1     = 5.
  l2     = 9.
  m1     = 13.
  m2     = 50.
  yi(1)  = 173.*rad
  yi(2)  = -86.*rad
  dydt(1) = -3.*rad
  dydt(2) = 35.*rad
  ms     = m1 + m2
  dtheta   = yi(1) - yi(2)
  yi(3)  = ms*l1*l1*dydt(1) + m2*l1*l2*dydt(2)*cos(dtheta)
  yi(4)  = m2*l2*l2*dydt(2) + m2*l1*l2*dydt(1)*cos(dtheta)
  ! Print initial conditions.
  write(1,*) t,yi(1),yi(2),yi(3),yi(4),l1*sin(yi(1)),-l1*cos(yi(1)),l1*sin(yi(1))+l2*sin(yi(2)),-l1*cos(yi(1))-l2*cos(yi(2))

  call cpu_time(start)
  do while (t < tmax)
    call rkck(d_pend_deriv,t,yi,yf,er,dt)
    ! Preparing next step.
    t = t + dt
    yi(:) = yf(:)
    ! Writing current results.
    write(1,*) t,yf(1),yf(2),yf(3),yf(4),l1*sin(yf(1)),-l1*cos(yf(1)),l1*sin(yf(1))+l2*sin(yf(2)),-l1*cos(yf(1))-l2*cos(yf(2))
  end do
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
end subroutine d_pend_stat

subroutine d_pend_adapt
  use d_pend
  use ode_int
  implicit none
  integer, parameter :: n = 4
  double precision, parameter :: rad = atan(1.)/45.
  double precision t, dt, dtmin, tmax, dtheta, der
  double precision y(n)
  double precision, dimension(size(y)) :: dydt, er
  double precision start, finish
  interface
    subroutine d_pend_deriv(t,y,dydt)
      implicit none
      double precision, intent(in) :: t, y(:)
      double precision, intent(out) :: dydt(:)
    end subroutine
  end interface

  open(unit=1, file = 'd_pend_adapt.plt')
  ! Initial conditions.
  t = 0.
  tmax = 50.
  dt     = 5.e-3
  dtmin  = 4.375e-3
  der    = 1.3125e-7
  l1     = 5.
  l2     = 9.
  m1     = 13.
  m2     = 50.
  y(1)  = 173.*rad
  y(2)  = -86.*rad
  dydt(1) = -3.*rad
  dydt(2) = 35.*rad
  ms     = m1 + m2
  dtheta   = y(1) - y(2)
  y(3)  = ms*l1*l1*dydt(1) + m2*l1*l2*dydt(2)*cos(dtheta)
  y(4)  = m2*l2*l2*dydt(2) + m2*l1*l2*dydt(1)*cos(dtheta)
  ! Print initial conditions.
  write(1,*) t,y(1),y(2),y(3),y(4),l1*sin(y(1)),-l1*cos(y(1)),l1*sin(y(1))+l2*sin(y(2)),-l1*cos(y(1))-l2*cos(y(2))

  call cpu_time(start)
  do while (t < tmax)
    call rkcka(d_pend_deriv,t,y,der,dt,dtmin)
    ! Writing current results.
    write(1,*) t,y(1),y(2),y(3),y(4),l1*sin(y(1)),-l1*cos(y(1)),l1*sin(y(1))+l2*sin(y(2)),-l1*cos(y(1))-l2*cos(y(2))
  end do
  call cpu_time(finish)
  print '("Time = ",f6.3," seconds.")',finish-start
end subroutine d_pend_adapt

subroutine d_pend_deriv(t,y,dydt)
! y(1) = theta1    dydt(1) = theta1'
! y(2) = theta2    dydt(2) = theta2'
! y(3) = p1        dydt(3) = p1'
! y(4) = p2        dydt(4) = p2'
use d_pend
implicit none
double precision  dtheta, den, c1, c2
double precision, intent(in) :: t, y(:)
double precision, intent(out):: dydt(:)
double precision, parameter :: g = 9.81
dtheta   = y(1) - y(2)
den   = m1 + m2*sin(dtheta)*sin(dtheta)
c1    = (y(3)*y(4)*sin(dtheta))/(l1*l2*den)
c2    = (l2*l2*m2*y(3)*y(3) + l1*l1*ms*y(4)*y(4) - l1*l2*m2*y(3)*y(4)*cos(dtheta))/(2*l1*l1*l2*l2*den*den)*sin(2*dtheta)
dydt(1) = (l2*y(3)-l1*y(4)*cos(dtheta))/(l1*l1*l2*den)
dydt(2) = (l1*ms*y(4)-l2*m2*y(3)*cos(dtheta))/(l1*l2*l2*m2*den)
dydt(3) = -ms*g*l1*sin(y(1)) - c1 + c2
dydt(4) = -m2*g*l2*sin(y(2)) + c1 - c2
end subroutine d_pend_deriv
