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
  type(Vec2D)  :: xi
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
  real(dp) :: coefficients(16)

  ! Initialise names.
  allocate(character(len=8) :: names2D(4))

  ! Test 2D.1: Initialise Fluid.
  ! Safety factor.
  safety = 1._dp
  ! Grid scale.
  scale = 1._dp
  ! 10x10 2D grid with sides of length 3 (each box is of size 1x1).
  whd2D = 5
  ! Density.
  rho = 0.1_dp
  ! Initial timestep.
  timestep = 0.005
  ! Names
  names2D(1)(:) = "Density."
  names2D(2)(:) = "Vx."
  names2D(3)(:) = "Vy."
  call FldSlv2D % Init(names2D, whd2D, rho, timestep, scale, safety)

  xi2D % x = [0.2_dp, .7_dp, 0.2_dp, 0.7_dp]
  initial_cond = [5.0_dp, 1.0_dp, 1.0_dp]
  xi % x = [0.5_dp, 0.5_dp]

  call FldSlv2D % InitFlow(xi2D, initial_cond)

  print*, "Slow Cubic interpolation: ", FldSlv2D % u(1) % Cerp(xi), new_line(new_line(""))


  coefficients = FldSlv2D % u(1) % BcC2D(xi)

  do i = 1, 16
    print*,i, coefficients(i)
  end do
  print*, new_line(new_line(""))

  print*, "Quick Cubic interpolation: ", FldSlv2D % u(1) % QuickCerp2D(xi, FldSlv2D % u(1) % BcC2D(xi)), new_line(new_line(""))
end program test_fluid
