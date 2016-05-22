program sfldy_main
  use sfldy
  real(dp) :: rho, h
  integer :: idx, i, j, k, sg, sx
  integer, parameter :: disp(24) = [0,1,0,1,0,1,0,1,0,0,1,1,&
                                    0,0,1,1,0,0,0,0,1,1,1,1]
  sx = 2
  sg = 2**sx
  idx = 0
  do i = 1, sg
    do j = 1, sx
      idx = idx + 1
      print*, "j loop", "index = ", i + 8*(j-1), "displacement = ", disp(i + 8*(j-1))
      if (mod(idx,sx) == 0) print*, "-------------"
      do k = 1, 3
        !if (j == k) cycle
        !print*, "k loop", i + 4*(k-1)
      end do
    end do
  end do


end program sfldy_main
