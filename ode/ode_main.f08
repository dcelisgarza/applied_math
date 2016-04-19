program ode_main
  use nrtype
  use derivatives
  use ode, only : basic_verlet, velocity_verlet
  implicit none
  real(dp) :: dvdt, t, xo(1), xn(1), h, vo, xov(2), xnv(2)

  h  = 0.15_dp
  xo = 0._dp
  t  = 0._dp
  vo = 15._dp
  open(unit=1, file='test1.dat')
  xn(1) = xo(1) + vo*h -0.5_dp*9.81_dp*h*h
  write(1,*) t, xo
  t = t + h
  write(1,*) t, xn
  do while (t < 3.0)
    call basic_verlet(constant_accel,t,xo,xn,h)
    t = t + h
    write(1,*) t, xn
  end do

  xov(1) = 0._dp
  xov(2) = 15._dp
  t      = 0._dp
  open(unit=1, file='test2.dat')
  write(1,*) t, xov(1), xov(2)
  do while (t < 3.0)
    call velocity_verlet(constant_accel,t,xov,xnv,h)
    t = t + h
    write(1,*) t, xnv(1), xnv(2)
    xov = xnv
  end do

end program ode_main
