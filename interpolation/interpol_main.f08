program interpol_main
  use numbers
  use interpolation
  real(dp) :: y, dy
  call poly_interpol([1._dp,2._dp],[2._dp,4._dp],1.5_dp,y,dy)
end program interpol_main
