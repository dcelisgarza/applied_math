program find_roots_main
  use numbers
  use continuous_functions
  use find_roots
  implicit none
  real(dp) :: root

  call bisection(polynomial,[1._dp,2._dp],1.d-6,20,root)
  print*, root

end program find_roots_main
