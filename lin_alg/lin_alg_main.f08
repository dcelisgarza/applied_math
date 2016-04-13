program lin_alg_main
  use lin_alg
  implicit none
  real(dp) :: a(4), b(3), c(3), x(4), y(4)

  a = [2._dp, 3._dp, 5._dp, 7._dp]
  c = [11._dp, 13._dp, 17._dp]
  b = [19._dp, 23._dp, 29._dp]
  y = [31._dp, 37._dp, 41._dp, 43._dp]

  call tridiag(c,a,b,y,x)
  print*, "using my solution"
  print*, x

  !print*, "--------------"



  call tridiag_nr(c,a,b,y,x)
  print*, "using nr's solution"
  print*, x
end program lin_alg_main
