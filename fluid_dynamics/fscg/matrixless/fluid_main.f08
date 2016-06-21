program fluid_main
  use fluid
  real(dp) :: celsiz

  real(dp), allocatable :: ost(:)
  integer, allocatable :: whd(:)
  type(FluidQty), allocatable :: FldQty(:)
  type(FluidSlv) :: FldSlv
  character(len=:), dimension(:), allocatable :: names

  integer :: i, j, k
  real(dp) :: l, dummy

  type(vec2d) :: xi

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

  ! Initialise FldSlv
  ! 2D
  deallocate(whd)
  allocate(whd(2))
  allocate(character(len = 8) :: names(size(whd) + 1))

  whd = 2 ! 2x2
  names(1)(1:len(names(1))) = "Pressure"
  names(2)(1:len(names(2))) = "Vel_X"
  names(3)(1:len(names(3))) = "Vel_Y"
  call FldSlv % init(names = names, whd = whd, scale = 10._dp)

  ! Testing pressure initialisation.
  print*, trim(names(1)(1:len(names(1))))
  print*, "name = ", FldSlv % rho % name
  print*, " # of cells = ", size(FldSlv % rho % old)
  print*, " value of old cells = ", FldSlv % rho % old
  print*, " # of cells = ", size(FldSlv % rho % new)
  print*, " value of new cells = ", FldSlv % rho % new
  print*, " whd = ", FldSlv % rho % whd
  print*, " ost = ", FldSlv % rho % ost
  print*, " size of cells = ", FldSlv % rho % celsiz

  print*, "-------------------------------"

  ! Testing Vel_X initialisation.
  print*, trim(names(2)(1:len(names(2))))
  print*, "name = ", FldSlv % u(1) % name
  print*, " # of cells = ", size(FldSlv % u(1) % old)
  print*, " value of old cells = ", FldSlv % u(1) % old
  print*, " # of cells = ", size(FldSlv % u(1) % new)
  print*, " value of new cells = ", FldSlv % u(1) % new
  print*, " whd = ", FldSlv % u(1) % whd
  print*, " ost = ", FldSlv % u(1) % ost
  print*, " size of cells = ", FldSlv % u(1) % celsiz

  print*, "-------------------------------"

  ! Testing Vel_Y initialisation.
  print*, trim(names(3)(1:len(names(3))))
  print*, "name = ", FldSlv % u(2) % name
  print*, " # of cells = ", size(FldSlv % u(2) % old)
  print*, " value of old cells = ", FldSlv % u(2) % old
  print*, " # of cells = ", size(FldSlv % u(2) % new)
  print*, " value of new cells = ", FldSlv % u(2) % new
  print*, " whd = ", FldSlv % u(2) % whd
  print*, " ost = ", FldSlv % u(2) % ost
  print*, " size of cells = ", FldSlv % u(2) % celsiz

  print*, "-------------------------------"

  ! Testing write and read.
  l = 0.
    do j = 1, nint(FldSlv % u(2) % whd(2))
      do i = 1, nint(FldSlv % u(2) % whd(1))
        l = l + 1.
        call FldSlv % u(2) % WriteVal(i,j,l)
      end do
    end do

  call FldSlv % u(2) % UpdateVals

    do j = 1, nint(FldSlv % u(2) % whd(2))
      do i = 1, nint(FldSlv % u(2) % whd(1))
        print*, " [i, j] = ", [i, j], " Value = ", FldSlv % u(2) % ReadVal(i,j)
      end do
    end do

    !xi % x = [0.5_dp,0.5_dp]
    print*, FldSlv % u(2) % Lerp(Vec2D([2._dp,0.5_dp]))

end program fluid_main
