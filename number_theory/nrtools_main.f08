program nrtools_main
  use nrtools
  implicit none
  integer :: n, k

  n = 4
  k = 5

  ! n!/(k!(n-k)!)

  print*, rcsv_bicoef(n,k)
end program nrtools_main
