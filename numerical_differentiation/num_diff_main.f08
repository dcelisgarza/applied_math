program num_diff_main
  use nrtype
  use num_diff
  real(dp) :: x, h
  real(dp) :: i
  real(dp) :: start, finish
  real(dp) :: aux

  !open(unit=1, file = 'test_numdif.dat')
  call cpu_time(start)
  do i = 0._dp, tau, 0.01_dp
    write(1,*) i, sin(i), cos(i),sifodi1(sine,i,0.0001_dp)
    , (cos(i)-sifodi1(sine,i,0.01_dp))*100._dp &
    ,eicedi1(sine,i,0.01_dp), (cos(i)-sifodi1(sine,i,0.01_dp))*100._dp
  end do
  call cpu_time(finish)
  print*, finish-start

contains
    function sine(x)
      real(dp), intent(in) :: x
      real(dp)             :: sine
      sine = sin(x)
    end function sine
end program num_diff_main
