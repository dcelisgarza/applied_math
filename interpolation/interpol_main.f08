program interpol_main
  use nrtype
  use interpolation
  implicit none
  real(dp) :: y, dy, x, grid(2,2,2), xyz(3)
  real(dp) :: xi(6), yi(6), yi2(6), yi3(6), yi4(6)
  integer :: i, j, idx(3)

  x  = -5._dp
  xi = [-5.0_dp,-2.2_dp,1._dp,2.2_dp,8.6_dp,9.3_dp]
  yi = [-8.7_dp,2._dp,4._dp,10.3_dp,-4.67_dp,8.6_dp]

  yi2 = rational(xi)

  call csplinec(xi,yi,yi3)
  call csplinec(xi,yi,yi4)

  open(unit = 1, file = "poly_int.dat")
  open(unit = 2, file = "rat_int.dat")
  open(unit = 3, file = "cspline_int.dat")
  open(unit = 4, file = "cspline_rat_int.dat")
  open(unit = 5, file = "trilinear.dat")

  do while (x < 9.3_dp)
    call pinex(xi,yi,x,y,dy)
    write(1,*) x, y

    call rinex(xi,yi2,x,y,dy)
    write(2,*) x, y, rational(x)

    call csplnei(xi, yi, yi3, x, y)
    write(3,*) x, y

    call csplnei(xi, yi2, yi4, x, y)
    write(4,*) x, y

    x = x + 0.1_dp
  end do
  close(1); close(2); close(3); close(4)

  ! Testing trilinear interpolation.
  grid = 1._dp
  xyz  = [1.5_dp,1.5_dp,1.5_dp]
  idx  = floor(xyz)
  ! Should print 1
  print*, cctrilinint(grid, xyz, idx, 1._dp)
  ! Test successful!

  grid(1,1,1) = 0.
  grid(1,2,1) = 0.
  grid(2,1,1) = 0.
  grid(1,1,2) = 0.
  grid(2,2,1) = 1.
  grid(2,1,2) = 1.
  grid(1,2,2) = 1.
  grid(2,2,2) = 1.
  xyz  = [1.9_dp,1.9_dp,1.9_dp]
  idx  = floor(xyz)
  ! Should print 0.972
  print*, cctrilinint(grid, xyz, idx, 1._dp)
  ! Test successful!

contains

  elemental function rational(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: rational

    rational = ( x*x*x - 2._dp * x ) / ( 2._dp*(x*x - 5._dp) )
  end function rational
end program interpol_main
