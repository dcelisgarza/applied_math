program networks_main
  use networks
  real(dp) :: ohd(6,6), dist(6)
  integer  :: prev(6)

  ohd = inf(1_i1)
  do i = 1, 6
    ohd(i,i) = 0
  end do

  ohd(1,2) = 40
  ohd(1,3) = 15
  ohd(2,3) = 20
  ohd(2,4) = 10
  ohd(2,5) = 25
  ohd(3,4) = 100
  ohd(2,6) = 6
  ohd(5,6) = 8

  ohd(2,1) = 40
  ohd(3,1) = 15
  ohd(3,2) = 20
  ohd(4,2) = 10
  ohd(5,2) = 25
  ohd(4,3) = 100
  ohd(6,2) = 6
  ohd(6,5) = 8

  call dijkstra(ohd, 1, dist, prev)
  print*, dist
  print*, prev
end program networks_main
