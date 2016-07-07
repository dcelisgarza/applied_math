program test_fluid
  use fluid
  implicit none
  ! Fluid solver class.
  type(FluidSlv) :: FldSlv2D, FldSlv3D
  ! Names of variables.
  character(len=:), dimension(:), allocatable :: names2D, names3D
  ! Dimensions of the simulation domains.
  integer :: whd2D(2), whd3D(3)
  ! Simulation parameters.
  real(dp) :: scale, safety
  ! Density.
  real(dp) :: rho, timestep
  ! Add in flow coordinates.
  type(Vec2D2) :: xi2D
  real(dp) :: initial_cond(3)
  ! Time
  real(dp) :: t
  ! Timers
  integer :: startvalues(8), endvalues(8), startcount, endcount
  real(4) :: countrate
  ! Counter
  integer :: i, j
  ! Arrays
  real(dp), allocatable :: u(:,:), prhs(:)
  real(dp) :: celsiz

  ! Initialise names.
  allocate(character(len=8) :: names2D(4))

  ! Test 2D.1: Initialise Fluid.
  ! Safety factor.
  safety = 1._dp
  ! Grid scale.
  scale = 1._dp
  ! 3x3 2D grid with sides of length 3 (each box is of size 1x1).
  whd2D = 3
  ! Density.
  rho = 0.1_dp
  ! Initial timestep.
  timestep = 0.005
  ! Names
  names2D(1)(:) = "Density."
  names2D(2)(:) = "Vx."
  names2D(3)(:) = "Vy."
  call FldSlv2D % Init(names2D, whd2D, rho, timestep, scale, safety)

  print*, " rho % name = ", FldSlv2D % rho % name, &
          new_line(" "), " Expected rho % name = ", names2D(1)(1:len(names2D(1)))
  print*, " u(1) % name = ", FldSlv2D % u(1) % name, &
          new_line(" "), " Expected u(1) % name = ", names2D(2)(1:len(names2D(2)))
  print*, " u(2) % name = ", FldSlv2D % u(2) % name, &
          new_line(" "), " Expected u(2) % name = ", names2D(3)(1:len(names2D(3)))

  print*, " Box Width = ", FldSlv2D % whd(1), &
          new_line(" "), " Expected width = ", whd2D(1)
  print*, " Box Height = ", FldSlv2D % whd(2), &
          new_line(" "), " Expected height = ", whd2D(2)

  print*, " Density = ", FldSlv2D % srho, &
          new_line(" "), " Expected Density = ", rho

  print*, " Cell size = ", FldSlv2D % celsiz, &
          new_line(" "), " Expected Cell size = ", scale/minval(whd2D)

  print*, " Safety factor = ", FldSlv2D % safety, &
          new_line(" "), " Expected Safety factor = ", safety

  print*, new_line(new_line(" ")), " FldSlv2D % prhs = ", FldSlv2D % prhs, &
          new_line(new_line(" ")), " Expected FldSlv2D % prhs = ", [(real(i)*0._dp, i = 1, product(whd2D))]

  print*, new_line(new_line(" ")), " FldSlv2D % p = ", FldSlv2D % prhs, &
          new_line(new_line(" ")), " Expected FldSlv2D % p = ", [(real(i)*0._dp, i = 1, product(whd2D))]

  print*, new_line(new_line(" ")), " FldSlv2D % u(1) % old = ", FldSlv2D % u(1) % old, &
          new_line(new_line(" ")), " Expected FldSlv2D % u(1) % old = ", [(real(i)*0._dp, i = 1, (whd2D(1) + 1)*whd2D(2))]

  print*, new_line(new_line(" ")), " FldSlv2D % u(1) % new = ", FldSlv2D % u(1) % new, &
          new_line(new_line(" ")), " Expected FldSlv2D % u(1) % new = ", [(real(i)*0._dp, i = 1, (whd2D(1) + 1)*whd2D(2))]

  print*, " Size(FldSlv2D % u(1) % old) = ", size(FldSlv2D % u(1) % old), &
          new_line(" "), " Expected Size(FldSlv2D % u(1) % old) = ", (whd2D(1) + 1)*whd2D(2)

  print*, " Size(FldSlv2D % u(1) % new) = ", size(FldSlv2D % u(1) % new), &
          new_line(" "), " Expected Size(FldSlv2D % u(1) % new) = ", (whd2D(1) + 1)*whd2D(2)

  print*, " FldSlv2D % u(1) % ost = ", FldSlv2D % u(1) % ost, &
          new_line(" "), " Expected FldSlv2D % u(1) % ost = ", [0._dp, 0.5_dp]

  print*, " FldSlv2D % u(1) % whd = ", FldSlv2D % u(1) % whd, &
          new_line(" "), " Expected FldSlv2D % u(1) % whd = ", [whd2D(1) + 1, whd2D(2)]

  print*, " FldSlv2D % u(1) % celsiz = ", FldSlv2D % u(1) % celsiz, &
          new_line(" "), " Expected FldSlv2D % u(1) % celsiz = ", scale/minval(whd2D)

  print*, new_line(new_line(" ")), " FldSlv2D % u(2) % old = ", FldSlv2D % u(2) % old, &
          new_line(new_line(" ")), " Expected FldSlv2D % u(2) % old = ", [(real(i)*0._dp, i = 1, (whd2D(2) + 1)*whd2D(1))]

  print*, new_line(new_line(" ")), " FldSlv2D % u(2) % new = ", FldSlv2D % u(2) % new, &
          new_line(new_line(" ")), " Expected FldSlv2D % u(2) % new = ", [(real(i)*0._dp, i = 1, (whd2D(2) + 1)*whd2D(1))]

  print*, " Size(FldSlv2D % u(2) % old) = ", size(FldSlv2D % u(2) % old), &
          new_line(" "), " Expected Size(FldSlv2D % u(2) % old) = ", (whd2D(2) + 1)*whd2D(1)

  print*, " Size(FldSlv2D % u(2) % new) = ", size(FldSlv2D % u(2) % new), &
          new_line(" "), " Expected Size(FldSlv2D % u(2) % new) = ", (whd2D(2) + 1)*whd2D(2)

  print*,  " FldSlv2D % u(2) % ost = ", FldSlv2D % u(2) % ost, &
          new_line(" "), " Expected FldSlv2D % u(2) % ost = ", [0.5_dp, 0._dp]

  print*, " FldSlv2D % u(2) % whd = ", FldSlv2D % u(2) % whd, &
          new_line(" "), " Expected FldSlv2D % u(2) % whd = ", [whd2D(1), whd2D(2) + 1]

  print*, " FldSlv2D % u(2) % celsiz = ", FldSlv2D % u(2) % celsiz, &
          new_line(" "), " Expected FldSlv2D % celsiz = ", scale/minval(whd2D)

  print*, new_line(new_line(" ")), " FldSlv2D % rho % new = ", FldSlv2D % rho % new, &
          new_line(new_line(" ")), " Expected FldSlv2D % rho % new = ", [(real(i)*0._dp, i = 1, product(whd2D))]

  print*, " Size(FldSlv2D % rho % old) = ", size(FldSlv2D % rho % old), &
          new_line(" "), " Expected Size(FldSlv2D % rho % old) = ", product(whd2D)

  print*, " Size(FldSlv2D % rho % new) = ", size(FldSlv2D % rho % new), &
          new_line(" "), " Expected Size(FldSlv2D % rho % new) = ", product(whd2D)

  print*,  " FldSlv2D % rho % ost = ", FldSlv2D % rho % ost, &
          new_line(" "), " Expected FldSlv2D % rho % ost = ", [0.5_dp, 0.5_dp]

  print*, " FldSlv2D % rho % whd = ", FldSlv2D % rho % whd, &
          new_line(" "), " Expected FldSlv2D % rho % whd = ", [whd2D(1), whd2D(2)]

  print*, " FldSlv2D % rho % celsiz = ", FldSlv2D % rho % celsiz, &
          new_line(" "), " Expected FldSlv2D % rho % celsiz = ", scale/minval(whd2D), new_line(new_line(" "))

  xi2D % x = [0._dp, 1._dp, 0._dp, 1._dp]
  initial_cond = [2.0_dp, 3.0_dp, 1._dp]

  CALL SYSTEM_CLOCK(startcount, countrate)
  do i = 1, 1000000
    call FldSlv2D % InitFlow(xi2D, initial_cond)
  end do
  CALL SYSTEM_CLOCK(endcount, countrate)
  print*, "InitFlow: It took me", endcount - startcount, "clicks (", (endcount - startcount)/countrate*1000," ms)"

  print*, new_line(new_line(" ")), "FldSlv2D % rho % old = ", new_line(" "), FldSlv2D % rho % old

  print*, new_line(new_line(" ")), "FldSlv2D % u(1) % old = ", new_line(" "), FldSlv2D % u(1) % old

  print*, new_line(new_line(" ")), "FldSlv2D % u(2) % old = ", new_line(" "), FldSlv2D % u(2) % old, new_line(new_line(" "))

  t = 0._dp
  CALL SYSTEM_CLOCK(startcount, countrate)
  do i = 1, 1000000
    call FldSlv2D % MaxTimeStep2D(t)
  end do
  CALL SYSTEM_CLOCK(endcount, countrate)
  print*, "MAXTIMESTEP: It took me", endcount - startcount, "clicks (", (endcount - startcount)/countrate*1000," ms)"
  print*, new_line(new_line(" ")), FldSlv2D % TimeStep


  CALL SYSTEM_CLOCK(startcount, countrate)
  do i = 1, 1000000
    call FldSlv2D % BuildRhs2D()
  end do
  CALL SYSTEM_CLOCK(endcount, countrate)
  print*, "BUILDRHS: It took me", endcount - startcount, "clicks (", (endcount - startcount)/countrate*1000," ms)"
  print*, new_line(new_line(" ")), " FldSlv2D % prhs = ", FldSlv2D % prhs

  CALL SYSTEM_CLOCK(startcount, countrate)
  call FldSlv2D % GSPSlv2D(1000000, 1d-5)
  CALL SYSTEM_CLOCK(endcount, countrate)
  !print*, endvalues - startvalues
  print*, "Gauss-Seidel: It took me", endcount - startcount, "clicks (", (endcount - startcount)/countrate*1000," ms)"

  call FldSlv2D % ApplyPressure2D()
  do j = 0, 2
    do i = 0, 2
      write(*,"(A,A,I2,A,I2,A,A,F5.3,A,F5.3,A)"), "(Vx(i,j), Vy(i,j) = ", "(", i,",", j, "), ", "(", FldSlv2D % u(1) % ReadVal(i,j)&
      ,", ", FldSlv2D % u(2) % ReadVal(i,j), ")"
    end do
  end do

  call FldSlv2D % rho % Advect2D(FluidQty2D([FldSlv2D % u]), FldSlv2D%timestep(1))
  call FldSlv2D % u(1) % Advect2D(FluidQty2D([FldSlv2D % u]), FldSlv2D%timestep(1))
  call FldSlv2D % u(2) % Advect2D(FluidQty2D([FldSlv2D % u]), FldSlv2D%timestep(1))
  do i = 0, 8
    print"(F12.10)", FldSlv2D % rho % new(i)
  end do
  do i = 0, 11
    print"(F12.10)", FldSlv2D % u(1) % new(i)
  end do
  do i = 0, 11
    print"(F12.10)", FldSlv2D % u(2) % new(i)
  end do

end program test_fluid
