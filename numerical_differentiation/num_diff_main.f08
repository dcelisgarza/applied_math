program num_diff_main
  use nrtype
  use num_diff
  real(dp) :: x, h
  real(dp) :: i
  real(dp) :: start, finish
  real(dp) :: aux
  integer, allocatable :: coefs(:)

  open(unit=1, file = 'test_numdif.dat')

  !do i = 0._dp, tau, 0.01_dp
  !  write(1,*) i, sin(i), cos(i),sifodi1(sine,i,0.0001_dp), (cos(i)-sifodi1(sine,i,0.01_dp))*100._dp &
  !  ,eicedi1(sine,i,0.01_dp), (cos(i)-sifodi1(sine,i,0.01_dp))*100._dp
  !end do
  !call cpu_time(finish)

  call cpu_time(start)

  allocate(coefs(2))
  coefs = dncoef1(1)

  print*, coefs

  do i = 0._dp, tau, 0.01_dp
    write(1,*) i, cos(i)-fdo1(sine, i, 0.01_dp, coefs)
  end do
  !print*, dncoef1(1), finish-start

contains
    function sine(x)
      real(dp), intent(in) :: x
      real(dp)             :: sine
      sine = sin(x)
    end function sine
end program num_diff_main
