program fluid_main
  use fluid
  implicit none
  integer  :: i, j, k, randint, bounds(6)
  real(dp) :: umax, tumax, r
  type(FluidCell), allocatable :: grid(:,:,:)
  integer :: grid_size(3)

  allocate( grid(3,3,3) )
  do i = 1, 3
    grid_size(i) = size(grid, dim = i)
  end do


!  call find_vmax(grid, umax)
!  print*, "SUBROUTINE find_vmax"
!  print*, "u max = ", umax

! bounds = [1,3,6,8,11,13]

!  doi: do i = 0, 4, 2
!    if (i < bounds(1) .or. i > bounds(2)) then
!      !print*, "i = ", i
!      cycle doi
!    end if
!    doj: do j = 5, 9, 2
!      if (j < bounds(3) .or. j > bounds(4)) then
!        !print*, "j = ", j
!        cycle doj
!      end if
!      dok: do k = 10, 14, 2
!        if (k < bounds(5) .or. k > bounds(6)) then
!          !print*, "k = ", k
!          cycle dok
!        end if
!        print*, "i = ", i, "j = ", j, "k = ", k
!      end do dok
!    end do doj
!  end do doi

! Successfully tested the part that updates the cells with liquid in them.

! No cells have markers.
! Expected result: cycles through all loops.
grid % marker = 0
call update_grid(1._dp, grid, grid_size)
! Test successful!

! Central grid has a liquid particle but is a solid.
! Expected result: Cycle through all loops where i, j, k /=  2.
! Note that the cell is a solid and do nothing regarding cell_type, and set cell_layer = -1.
grid(2,2,2) % marker = 1
grid(2,2,2) % cell_type = 3
call update_grid(1._dp, grid, grid_size)
print*, 1, grid(2,2,2) % cell_type, grid(2,2,2) % layer
! Test successful!

! Central grid has a liquid particle and is solid (i = 1 -> gas, i = 2 -> liquid).
! Expected result: cell_type = 2 (liquid), layer = 0.
grid(2,2,2) % marker = 1
do i = 1, 2
  grid(2,2,2) % cell_type = i
  call update_grid(1._dp, grid, grid_size)
  print*, i, grid(2,2,2) % cell_type, grid(2,2,2) % layer
end do
! Test successful!

! Test the buffer zone.

! Central cell is i, buffer cells are j, buffer cells are NOT initial.
! Expected result: buffers will be air, buffer cells layer /= -1.
doi: do i = 1, 2
  doj: do j = 1, 2
    grid % cell_type = 3
    grid % marker = 0
    grid(2,2,2) % marker = i
    grid(2,2,2) % cell_type = 2
    grid(2,2,2) % buffer % initial = .false.
    grid(2,2,2) % buffer % cell_type = j
    call update_grid(1._dp, grid, grid_size)
    print*, "i = ", i, "j = ", j
    print*, "buffer % cell_type ="
    print*, grid(2,2,2) % buffer % cell_type
    print*, "buffer % layer ="
    print*, grid(2,2,2) % buffer % layer
  end do doj
end do doi
! Test successful!

! Central cell is i, buffer cells are j, buffer cells are initial. Buffer cells are within bounds.
! Expected result: buffers will be air, buffer cells layer /= -1.
doi2: do i = 1, 2
  doj2: do j = 1, 2
    grid % cell_type = 3
    grid % marker = 0
    grid(2,2,2) % marker = i
    grid(2,2,2) % cell_type = 2
    grid(2,2,2) % buffer % initial = .true.
    grid(2,2,2) % buffer % cell_type = j
    call update_grid(1._dp, grid, grid_size)
    print*, "i = ", i, "j = ", j
    print*, "buffer % cell_type ="
    print*, grid(2,2,2) % buffer % cell_type
    print*, "buffer % layer ="
    print*, grid(2,2,2) % buffer % layer
    print*, "buffer % initial ="
    print*, grid(2,2,2) % buffer % initial
  end do doj2
end do doi2
! Test successful!

! Central cell is the left-top-front grid(1,3,3). Therefore the buffers to the left grid(0,3,3), above grid(1,4,3) and in front grid(1,3,4) are out of bounds and are therefore solid.
print*, "TESTING BOUNDS"
doi3: do i = 1, 2
  doj3: do j = 1, 2
    grid % cell_type = 3
    grid % marker = 0
    grid(1,3,3) % marker = i
    grid(1,3,3) % cell_type = 2
    grid(1,3,3) % buffer % initial = .true.
    grid(1,3,3) % buffer % cell_type = j
    call update_grid(1._dp, grid, grid_size)
    print*, "i = ", i, "j = ", j
    print*, "buffer % cell_type ="
    print*, grid(1,3,3) % buffer % cell_type
    print*, "buffer % layer ="
    print*, grid(1,3,3) % buffer % layer
    print*, "buffer % initial ="
    print*, grid(1,3,3) % buffer % initial
  end do doj3
end do doi3

! All cells are solid and have layer = -1, all buffers are not initial.
! Expected result: No buffer zone is created, all buffer zones are initialised (initial = .true.) i.e. "deleted"
do i = 1, 6
  grid % buffer(i) % initial = .false.
end do
grid % cell_type = 3
call update_grid(1._dp, grid, grid_size)
do i = 1, 6
  print*, grid % buffer(i) % initial
end do
! Test successful.

end program fluid_main
