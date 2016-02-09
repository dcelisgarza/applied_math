program interpol_main
  use nrtype
  use interpolation
  real(dp) :: y, dy, x
  real(dp) :: xi(6), yi(6)
  x  = -5._dp
  xi = [-5.0_dp,-2.5_dp,1._dp,2._dp,8.6_dp,9.3_dp]
  yi = [-8.7_dp,2._dp,4._dp,10.3_dp,-4.67_dp,8.6_dp]
  open(unit=1,file='poly_int.dat')
  do while (x < 9.3_dp)
    call poly_interpol(xi,yi,x,y,dy)
  write(1,*) x, y
  x = x + 0.1_dp
  end do
end program interpol_main
