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

  allocate(coefs(3))
  coefs = dncoef1(2)

  do i = 0._dp, tau, 0.01_dp
    write(1,*) i, fdo1(sine, i, 0.01_dp, coefs), bdo1(sine, i, 0.01_dp, coefs), cdo1(sine, i, 0.01_dp, coefs)
  end do

  !deallocate(coefs)
  !aux = 0.
  !allocate(coefs(3))
  !coefs = dncoefn(1,2)
  !print*, aux,aux,aux,aux, dncoef1(1)
  !print*, aux,aux,aux,-dncoef1(2)/2
  !print*, aux,aux, dncoef1(3)/3
  !print*, aux,-dncoef1(4)/4
  !print*, dncoef1(5)/5
  !print*, "--------"
!  print*, aux,aux, aux, dncoef1(2)
!  print*, aux,aux, -dncoef1(3)
!  print*, aux, dncoef1(4)
!  print*, dncoef1(5)
!  print*, dncoef1(6)
  !print*, dncoefn(2,3)
  print*, dncoefn(2, 4, [0._dp,1._dp,2._dp,3._dp,4._dp], 0._dp)


contains
    function sine(x)
      real(dp), intent(in) :: x
      real(dp)             :: sine
      sine = sin(x)
    end function sine
end program num_diff_main
