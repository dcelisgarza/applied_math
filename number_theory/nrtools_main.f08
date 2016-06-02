program nrtools_main
  use nrtools
  implicit none
  integer :: n, k
  real(sp):: i
  real(sp):: a, b
  complex(dpc) :: test
  integer :: j
  real(dp) :: start, finish
  real(dp) :: rln_gamma

  test = 0.


  !n = 4
  !k = 5

  ! n!/(k!(n-k)!)

  !print*, rcsv_bicoef(n,k)

  !test = (1._dp,0._dp)

  !open(unit=1, file='gammatest.dat')
  !do i = 1., 5., 0.1
    !write(1,*) i-1, rcsv_gamma(test)
  !  test = test + 0.1
  !end do

  test = (142._dp, 0._dp)

  call cpu_time(start)
  do j = 1, 1000000
    rln_gamma = real(rcsv_ln_gamma(test))
  end do
  call cpu_time(finish)

  print*, "rcsv_ln_gamma execution time = ", finish - start

  call cpu_time(start)
  do j = 1, 1000000
    rln_gamma = log(real(rcsv_gamma(test)))
  end do
  call cpu_time(finish)

  print*, "log(rcsv_gamma) execution time = ", finish - start


end program nrtools_main
