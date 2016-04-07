program interpol_main
  use nrtype
  use interpolation
  implicit none
  real(dp) :: y, dy, x
  real(dp) :: xi(6), yi(6), yi2(6)
  integer :: i

  x  = -5._dp
  xi = [-5.0_dp,-2.2_dp,1._dp,2.2_dp,8.6_dp,9.3_dp]
  yi = [-8.7_dp,2._dp,4._dp,10.3_dp,-4.67_dp,8.6_dp]

  yi2 = rational(xi)

  open(unit = 1, file = "poly_int.dat")
  open(unit = 2, file = "rat_int.dat")
  do while (x < 9.3_dp)
    call pinex(xi,yi,x,y,dy)
    write(1,*) x, y

    call rinex(xi,yi2,x,y,dy)
    write(2,*) x, y, rational(x)

    x = x + 0.1_dp
  end do
  close(1); close(2)




contains
  elemental function rational(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: rational

    rational = ( x*x*x - 2._dp * x ) / ( 2._dp*(x*x - 5._dp) )
  end function rational
end program interpol_main
