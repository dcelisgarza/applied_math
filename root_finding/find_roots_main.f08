program find_roots_main
  use nrtype
  use continuous_functions
  use find_roots
  implicit none
  real(dp) :: root
  real(dp) :: interval(2)
  logical  :: success
  real(dp) :: interval_seg(2,3)
  integer  :: n_root

  interval = [-11_dp,6_dp]
  n_root = 3
  call inbracket(polynomial,interval,interval_seg,6,n_root)
  print*, interval_seg, n_root

  !call outbracket(polynomial,interval,success,20,1.5_dp)

  !interval = [-2_dp,5_dp]
  !call bisection(polynomial,interval,1.d-8,30,root)
  !print*, success, interval, root

end program find_roots_main
