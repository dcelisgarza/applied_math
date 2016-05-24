program fluid_main
  use fluid
  implicit none
  integer  :: i, j, k, randint, bounds(6)
  real(dp) :: umax, tumax, r
  type(FluidCell), allocatable :: grid(:,:,:)

  allocate( grid(3,2,1) )

  call find_vmax(grid, umax)
  print*, "SUBROUTINE find_vmax"
  print*, "u max = ", umax

 bounds = [1,3,6,8,11,13]

  doi: do i = 0, 4, 2
    if (i < bounds(1) .or. i > bounds(2)) then
      !print*, "i = ", i
      cycle doi
    end if
    doj: do j = 5, 9, 2
      if (j < bounds(3) .or. j > bounds(4)) then
        !print*, "j = ", j
        cycle doj
      end if
      dok: do k = 10, 14, 2
        if (k < bounds(5) .or. k > bounds(6)) then
          !print*, "k = ", k
          cycle dok
        end if
        print*, "i = ", i, "j = ", j, "k = ", k
      end do dok
    end do doj
  end do doi
end program fluid_main
