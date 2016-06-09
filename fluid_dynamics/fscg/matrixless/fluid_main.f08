program fluid_main
  use fluid
  real(dp) :: celsiz

  real(dp), allocatable :: ost(:)
  integer, allocatable :: whd(:)
  type(FluidQty), allocatable :: FldQty(:)

  integer :: i, j, k
  real(dp) :: l, dummy

  allocate(FldQty(2))

  ! 2D
  allocate(ost(2), whd(2))

  ! Testing initialisation of the type.
  whd = 2 ! 2x2 grid
  ost = 1.
  celsiz = 1.
  call FldQty(1) % init("BABYMETAL", whd, ost, celsiz)
  print*, " name = ", FldQty(1) % name
  print*, " # of cells = ", size(FldQty(1) % old)
  print*, " value of old cells = ", FldQty(1) % old
  print*, " # of cells = ", size(FldQty(1) % new)
  print*, " value of new cells = ", FldQty(1) % new
  print*, " whd = ", FldQty(1) % whd
  print*, " ost = ", FldQty(1) % ost
  print*, " size of cells = ", FldQty(1) % celsiz

  ! Testing write and read.

  l = 0.
  do j = 1, size(whd)
    do i = 1, size(whd)
      l = l + 1.
      call FldQty(1) % WriteVal(i,j,l)
    end do
  end do

  call FldQty(1) % UpdateVals

  do j = 1, size(whd)
    do i = 1, size(whd)
      print*, " [i, j] = ", [i, j], " Value = ", FldQty(1) % ReadVal(i,j)
    end do
  end do

  ! 3D
  deallocate(ost, whd)
  allocate(ost(3), whd(3))

  ! Testing initialisation of the type.
  whd = 3 ! 3x3x3 grid
  ost = 1.
  celsiz = 1.
  call FldQty(2) % init("BABYMETAL", whd, ost, celsiz)
  print*, " name = ", FldQty(2) % name
  print*, " # of cells = ", size(FldQty(2) % old)
  print*, " value of old cells = ", FldQty(2) % old
  print*, " # of cells = ", size(FldQty(2) % new)
  print*, " value of new cells = ", FldQty(2) % new
  print*, " whd = ", FldQty(2) % whd
  print*, " ost = ", FldQty(2) % ost
  print*, " size of cells = ", FldQty(2) % celsiz

  ! Testing write and read.

  l = 0.
  do k = 1, size(whd)
    do j = 1, size(whd)
      do i = 1, size(whd)
        l = l + 1.
        !print*, l
        call FldQty(2) % WriteVal(i,j,k,l)
      end do
    end do
  end do

  call FldQty(2) % UpdateVals

  do k = 1, size(whd)
    do j = 1, size(whd)
      do i = 1, size(whd)
        print*, " [i, j, k] = ", [i, j, k], " Value = ", FldQty(2) % ReadVal(i,j,k)
      end do
    end do
  end do

end program fluid_main
