! Generalised Fortran version of original c++ code by TunnaBrain (Benedikt Bitterli) https://github.com/tunabrain/incremental-fluids.

module fluid
  use nrtype

  type FluidQty
    real(dp)        , allocatable :: old(:) ! Old value.
    real(dp)        , allocatable :: new(:) ! New value.
    real(dp)        , allocatable :: ost(:) ! Offset.
    real(dp)                      :: celsiz ! Cell size (area or volume).
    integer         , allocatable :: whd(:) ! Width, height, depth of the simulation domain.
    character(len=:), allocatable :: name   ! Variable name.
    !dir$ attributes align:8 :: old
    !dir$ attributes align:8 :: new
    !dir$ attributes align:8 :: ost
    !dir$ attributes align:4 :: whd
  contains
    ! Initialise type.
    procedure :: Init => InitFluidQty
    ! Add fluid flow.
    procedure :: InitFlow2D
    procedure :: InitFlow3D
    generic   :: InitFlow => InitFlow2D, InitFlow3D
    ! Read values.
    procedure :: ReadVal2D
    procedure :: ReadVal3D
    generic   :: ReadVal => ReadVal2D, ReadVal3D
    ! Write values.
    procedure :: WriteNewVal2D
    procedure :: WriteNewVal3D
    generic   :: WriteNewVal => WriteNewVal2D, WriteNewVal3D
    procedure :: WriteOldVal2D
    procedure :: WriteOldVal3D
    generic   :: WriteOldVal => WriteOldVal2D, WriteOldVal3D
    ! Update values.
    procedure :: UpdateVals
    ! Linear interpolate.
    procedure :: Lerp1D
    procedure :: Lerp2D
    procedure :: Lerp3D
    generic   :: Lerp => Lerp1D, Lerp2D, Lerp3D
    ! Integrate.
    procedure :: FEulerVel2D
    procedure :: FEulerVel3D
    generic   :: FEuler => FEulerVel2D, FEulerVel3D
    ! Advect.
    procedure :: Advect2D
    procedure :: Advect3D
    generic   :: Advect => Advect2D, Advect3D
  end type FluidQty

  type FluidSlv
    type(FluidQty) , allocatable :: u(:)     ! Velocities
    type(FluidQty)               :: rho      ! Density array
    real(dp)       , allocatable :: prhs(:)  ! Right hand side of pressure equation.
    real(dp)       , allocatable :: p(:)     ! Pressure solution.
    real(dp)                     :: TimeStep(3) ! Timesteps
    real(dp)                     :: celsiz   ! Cell size.
    real(dp)                     :: srho     ! Scalar (fluid) density.
    real(dp)                     :: safety   ! Safety factor for timestep calculation
    integer        , allocatable :: whd(:)   ! Width, height, depth
    !dir$ attributes align:8 :: u
    !dir$ attributes align:8 :: prhs
    !dir$ attributes align:8 :: p
    !dir$ attributes align:4 :: whd
  contains
    ! Initialise type.
    procedure :: Init => InitFluid
    ! Add flow.
    procedure :: InitFlowQties2D
    procedure :: InitFlowQties3D
    generic   :: InitFlow => InitFlowQties2D, InitFlowQties3D
    ! Max time step.
    procedure :: MaxTimestep2D
    procedure :: MaxTimestep3D
    ! Build right hand side of pressure solution.
    procedure :: BuildRhs2D
    procedure :: BuildRhs3D
    ! Solve pressure.
    procedure :: GSPSlv2D
    procedure :: GSPSlv3D
    ! Apply pressure.
    procedure :: ApplyPressure2D
    procedure :: ApplyPressure3D
    ! Write values.
    procedure :: WriteArray2D
    procedure :: WriteArray3D
    generic   :: WriteArray => WriteArray2D, WriteArray3D
    !Read values.
    procedure :: ReadArray2D
    procedure :: ReadArray3D
    generic   :: ReadArray => ReadArray2D, ReadArray3D
    ! Write to file.
    procedure :: WriteToFile
  end type FluidSlv

  type FluidQty3D
    type(FluidQty) :: FldQty(3)
  end type FluidQty3D

  type FluidQty2D
    type(FluidQty) :: FldQty(2)
  end type FluidQty2D

  type Vec3D2
    real(dp) :: x(6)
  end type Vec3D2

  type Vec2D2
    real(dp) :: x(4)
  end type Vec2D2

  type Vec3D
    real(dp) :: x(3)
  end type Vec3D

  type Vec2D
    real(dp) :: x(2)
  end type Vec2D

contains

!==============================================================================!
! Tested.
  subroutine InitFluid(FldSlv, names, whd, rho, timestep, scale, safety)
    implicit none
    class(FluidSlv) , intent(out)  :: FldSlv ! Fluid
    integer         , intent(in)   :: whd(:) ! Width (i), height (j), depth (k)
    real(dp), dimension(size(whd)) :: ost    ! Offset
    real(dp), dimension(size(whd)) :: TMPost
    real(dp)        , intent(in)   :: rho
    real(dp)        , intent(in)   :: timestep
    real(dp), optional, intent(in) :: scale
    real(dp), optional, intent(in) :: safety
    integer, dimension(size(whd))  :: TMncells ! Number of cells for each different quantity
    integer                        :: ndim   ! Number of dimensions
    integer                        :: ncells ! Number of cells
    integer                        :: i      ! Counter
    character(len=*), dimension(:), intent(in) :: names
    !DIR$ ASSUME_ALIGNED ost: 8
    !DIR$ ASSUME_ALIGNED TMPost: 8
    !DIR$ ASSUME_ALIGNED whd: 4
    !DIR$ ASSUME_ALIGNED TMncells: 4

    ! Number of dimensions.
    ndim   = size(whd)
    ncells = product(whd) - 1
    !dir$ vector aligned
    ost    = 0.5

    ! Cell size.
    if(present(scale) .eqv. .true.) then
      FldSlv % celsiz = scale/minval(whd)
    else
      FldSlv % celsiz = 1./minval(whd)
    endif

    ! Initialise velocities.
    ! Staggered on the grid.
    allocate(FldSlv % u(ndim))
    iv: do i = 1, ndim
      !dir$ vector aligned
      TMncells    = whd
      TMncells(i) = whd(i) + 1
      !dir$ vector aligned
      TMPost    = ost
      TMPost(i) = ost(i) - 0.5
      call FldSlv % u(i) % init(trim(names(i+1)(1:len(names(i+1)))), TMncells, TMPost, FldSlv % celsiz)
    end do iv

    ! Initialise density.
    ! Centered in the grid.
    call FldSlv % rho % init(trim(names(1)(1:len(names(1)))), whd, ost, FldSlv % celsiz)

    ! Allocate the right hand side of pressure solve.
    allocate(FldSlv % prhs(0:ncells))
    !dir$ loop count min(512)
    !dir$ vector aligned
    FldSlv % prhs = 0._dp

    ! Allocate the pressure solution and give initial guess.
    allocate(FldSlv % p(0:ncells))
    !dir$ loop count min(512)
    !dir$ vector aligned
    FldSlv % p = 0._dp

    ! Set initial time step.
    FldSlv % TimeStep(1) = timestep
    if(present(safety) .eqv. .true.) then
      FldSlv % TimeStep(2) = timestep * safety * 1d-3
      FldSlv % TimeStep(3) = timestep * safety * 1d+3
    else
      FldSlv % TimeStep(2) = timestep * 0.8_dp * 1d-3
      FldSlv % TimeStep(3) = timestep * 0.8_dp * 1d+3
    end if

    ! Set the fluid density.
    !dir$ vector aligned
    FldSlv % srho = rho

    ! Safety factor.
    if(present(safety) .eqv. .true.) then
      FldSlv % safety = safety
    else
      FldSlv % safety = 0.8_dp
    endif

    ! Allocate whd (simulation bounds).
    allocate (FldSlv % whd(ndim))
    !dir$ vector aligned
    FldSlv % whd = whd
  end subroutine InitFluid
!==============================================================================!

!==============================================================================!
! Tested
  subroutine InitFluidQty(FldQty, name, whd, ost, celsiz)
    implicit none
    class(FluidQty),  intent(out) :: FldQty ! Fluid Quantity
    real(dp),         intent(in)  :: ost(:) ! Offset
    real(dp),         intent(in)  :: celsiz ! Cell size (area or volume)
    integer,          intent(in)  :: whd(:) ! Width (i), height (j), depth (k)
    integer                       :: ncells ! # of cells = prod(whd)
    character(len=*), intent(in)  :: name   ! Variable name.
    !DIR$ ASSUME_ALIGNED ost: 8
    !DIR$ ASSUME_ALIGNED whd: 4

    ! Check that size(whd) == size(ost).
    check_size: if ( size(whd) /= size(ost) ) then
      write(*,*) " Error: fluid.f08: InitFluidQty: check_size: Variable name = ", name, "."
      write(*,*) " Size(whd) == size(ost); size(whd) = ", size(whd), " size(ost) = ", size(ost)
      stop
    end if check_size

    ! Allocate values.
    ncells = product(whd) - 1
    allocate(FldQty % old(0:ncells), FldQty % new(0:ncells))
    !dir$ loop count min(512)
    !dir$ vector aligned
    FldQty % old = 0._dp
    !dir$ loop count min(512)
    !dir$ vector aligned
    FldQty % new = 0._dp

    ! Offset
    allocate(FldQty % ost(size(ost)))
    !dir$ vector aligned
    FldQty % ost = ost

    ! Cell size
    FldQty % celsiz = celsiz

    ! Width, height, depth
    allocate(FldQty % whd(size(whd)))
    !dir$ vector aligned
    FldQty % whd = whd

    ! Allocate name.
    allocate(character(len=len(name)):: FldQty % name)
    FldQty % name = name

  end subroutine InitFluidQty
!==============================================================================!

!==============================================================================!
  ! Initialise flow.
  ! Adds the respective quantities inside the box bounded by the coordinates xi%x.
  ! Tested: 2D
  subroutine InitFlowQties2D(FldSlv, xi, values)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    type(Vec2D2)   , intent(in)    :: xi        ! Coordinates for the vertices of the 2D box.
    real(dp)       , intent(in)    :: values(3) ! value(1) := density, value(2) = Vx, value(3) = Vy.

    call FldSlv % u(1:2) % InitFlow(xi, values(1:2))
    call FldSlv % rho    % InitFlow(xi, values(3)  )
  end subroutine InitFlowQties2D

  pure subroutine InitFlowQties3D(FldSlv, xi, values)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    type(Vec3D2)   , intent(in)    :: xi        ! Coordinates for the vertices of the 3D box.
    real(dp)       , intent(in)    :: values(4) ! value(1) := density, value(2) = Vx, value(3) = Vy, value(4) = Vz

    call FldSlv % u(1:3) % InitFlow(xi, values(1:3))
    call FldSlv % rho    % InitFlow(xi, values(4)  )
  end subroutine InitFlowQties3D
!==============================================================================!

!==============================================================================!
  ! Add relevant quantities to the fluid simulation.
  ! Tested: 2D.
  elemental subroutine InitFlow2D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec2D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1
    real(dp)       , intent(in)    :: value
    integer                        :: ixi(4), x, y, AuxY, idx

    ! (x0,y1) ----------- (x1,y1)
    !         |         |
    !         |         |
    !         |         |
    ! (x0,y0) ----------- (x1,y0)

    ! Set up the simulation bounds.
    !dir$ vector aligned
    bounds: do x = 1, 2
      y = 2*x
      !dir$ vector aligned
      ixi(y-1: y) = aint( xi % x(y-1: y)  / FldQty % celsiz - FldQty % ost(x) )
    end do bounds

    !dir$ vector aligned
    do y = min(ixi(4), FldQty % whd(2)) - 1, max(ixi(3), 0), -1
      ! Loop through x.
      ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
      AuxY = FldQty % whd(1)*y
      !dir$ loop count min(512)
      !dir$ vector aligned
      FldQty % old(min(ixi(2), FldQty % whd(1)) - 1 + AuxY: max(ixi(1), 0) + AuxY: -1) = value
    end do

  ! [(AuxY*FldQty % whd(1), AuxY = min(ixi(4), FldQty % whd(2)) - 1, max(ixi(3), 0), -1)]

  end subroutine InitFlow2D
!==============================================================================!

!==============================================================================!
  elemental subroutine InitFlow3D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec3D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1, xi%x(5) = z0, xi%x(6) = z1
    real(dp)       , intent(in)    :: value
    integer                        :: ixi(6), x, y, z, AuxZ, AuxY, idx

    ! This is the simulation domain.
    !
    !    (x0,y1,z1)   ------------- (x1,y1,z1)
    !                /|          /|
    !               / |         / |
    !              /  |        /  |
    ! (x0,y1,z0)  -------------   | (x1,y1,z0)
    !             |   |       |   |
    !  (x0,y0,z1) |   --------|--- (x1,y0,z1)
    !             |  /        |  /
    !             | /         | /
    !             |/          |/
    ! (x0,y0,z0)  ------------- (x1,y0,z0)

    bounds: do x = 1, 3
      y = 2*x
      !dir$ vector aligned
      ixi(y-1 : y) = aint( (xi % x(y-1 : y) ) / FldQty % celsiz - FldQty % ost(x) )
    end do bounds

    ! Loop through z.
    ! Start from the maximum value between the desired z0 and 1 (for the cell index) up to the minimum value between the desired z1 and the simulation depth.
!    lz: do z = min(ixi(6), FldQty % whd(3)) - 1, max(ixi(5), 0), -1
      ! Loop through y.
      ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
!      AuxZ = FldQty % whd(1)*FldQty % whd(2)*z
!      ly: do y = min(ixi(4), FldQty % whd(2)) - 1, max(ixi(3), 0), -1
        ! Loop through x.
        ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
!        AuxY = FldQty % whd(1)*y
!        lx: do x = min(ixi(2), FldQty % whd(1)) - 1, max(ixi(1), 0), -1
          ! Write value to old value.
!          idx = x + AuxY + AuxZ
          !if ( abs(FldQty % old(idx)) < abs(value) ) FldQty % old(idx) = value
!          FldQty % old(idx) = value
!        end do lx
!      end do ly
!    end do lz

    !dir$ vector aligned
    do z = min(ixi(6), FldQty % whd(3)) - 1, max(ixi(5), 0), -1
      ! Loop through y.
      ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
      AuxZ = FldQty % whd(1)*FldQty % whd(2)*z
      do y = min(ixi(4), FldQty % whd(2)) - 1, max(ixi(3), 0), -1
        ! Loop through x.
        ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
        AuxY = FldQty % whd(1)*y
        !dir$ loop count min(512)
        !dir$ vector aligned
        FldQty % old(min(ixi(2), FldQty % whd(1)) - 1 + AuxY + AuxZ: max(ixi(1), 0) + AuxY + AuxZ: -1) = value
      end do
    end do
  end subroutine InitFlow3D
!==============================================================================!

!==============================================================================!
! Tested: 2D.
  subroutine BuildRhs(FldSlv)
    ! Build Rhs wrapper subroutine.
    ! The actual subroutines build the pressure solution via the negative of the divergence as a forward finite difference (truly horrible). This is not the best way to go about doing this, but it will be fixed in the future.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    integer, allocatable, save     :: ndim ! Number of dimensions.

    ! Allocate ndim and set its value to the number of dimensions.
    if (.not. allocated(ndim)) allocate(ndim, source = size(FldSlv % whd))

    ! Choose 3D or 2D.
    C2D3D: if (ndim == 2) then
      call FldSlv % BuildRhs2D
    else C2D3D
      call FldSlv % BuildRhs3D
    end if C2D3D
  end subroutine BuildRhs
  !==============================================================================!

  !==============================================================================!
  ! Tested: 2D
  pure subroutine BuildRhs2D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale     ! Scaling factor.
    integer                        :: idx, y, x ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz
    idx = 0
    ly: do y = 0, FldSlv % whd(2) - 1
      lx: do x = 0, FldSlv % whd(1) - 1
        FldSlv % prhs(idx) = -scale * (                                         &
                                             FldSlv % u(1) % ReadVal(x+1, y  )  &
                                       +     FldSlv % u(2) % ReadVal(x  , y+1)  &
                                       - sum(FldSlv % u    % ReadVal(x , y  ))  &
                                      )
      idx = idx + 1
      end do lx
    end do ly
  end subroutine BuildRhs2D
  !==============================================================================!

  !==============================================================================!
  pure subroutine BuildRhs3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale
    integer                        :: idx, z, y, x ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz
    idx = 0
    lz: do z = 0, FldSlv % whd(3) - 1
      ly: do y = 0, FldSlv % whd(2) - 1
        lx: do x = 0, FldSlv % whd(1) - 1
          FldSlv % prhs(idx) = -scale * (                                              &
                                               FldSlv % u(1) % ReadVal(x+1, y  , z  )  &
                                         +     FldSlv % u(2) % ReadVal(x  , y+1, z  )  &
                                         +     FldSlv % u(3) % ReadVal(x  , y  , z+1)  &
                                         - sum(FldSlv % u    % ReadVal(x  , y  , z  )) &
                                        )
          idx = idx + 1
        end do lx
      end do ly
    end do lz
  end subroutine BuildRhs3D
  !==============================================================================!

  !==============================================================================!
  ! Tested: 2D.
  pure subroutine GSPSlv2D(FldSlv, MaxIt, DErr)
    ! Gauss-Seidel pressure solve.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: DErr   ! Desired error threshold.
    integer        , intent(in)    :: MaxIt  ! Maximum number of iterations.
    real(dp)                       :: scale, Err, err2 ! Scaling factor, iteration error.
    real(dp)                       :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                       :: NewP ! New pressure value.
    integer                        :: AuxWHD(2)
    integer                        :: iter, idx, x, y ! Iteration, counters, index
    !dir$ attributes align:8 :: AuxWHD

    ! Set scale.
    scale  = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
    !dir$ vector aligned
    AuxWHD = FldSlv % whd - 1

    ! ITERation Loop.

    iterl: do iter = 1, MaxIt
      ! Loop over y.
      Err = 0._dp
      idx = 0
      ly: do y = 0, AuxWHD(2)
        ! Loop over x.
        lx: do x = 0, AuxWHD(1)
          Diag    = 0._dp
          OffDiag = 0._dp

          ! Building the matrix as a five-point stencil. Solid borders, i.e. no fluid is outside of the simulation bounds.

          ! The pressure contribution for each cell depends on its 4 nearest orthogonal neighbours.
          ! For example, cell a22 will receive contributions from cells a12, a23, a33 and a21.
          ! Boundary cells will only receive contributions from the orthogonally placed cells that are found within the simulation bounds. For example, cell a13 will only receive contributions from cells a12 and a23, while cell a21 from cells a11, a22 and a31. This leads to the series of if statements found herein.

          ! y-axis
          ! ^
          ! |__> x-axis
          !
          ! -------------------
          ! | a13 | a23 | a33 |
          ! -------------------
          ! | a12 | a22 | a32 |     --------->   [a11,a21,a31,a12,a22,a32,a13,a23,a33]
          ! -------------------
          ! | a11 | a21 | a31 |
          ! -------------------scale

          ! Exclude Lower X Bound (Exclude cells a11, a12, a13 in the diagram).
          elxb: if (x > 0) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx - 1)
          end if elxb

          ! Exclude Lower X Bound (Exclude cells a13, a23, a33 in the diagram).
          euxb: if (x < AuxWHD(1)) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx + 1)
          end if euxb

          ! Exclude Upper Y Bound (Exclude cells a11, a12, 9 in the diagram).
          elyb: if (y > 0) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1))
          end if elyb

          ! Exclude Lower Y Bound (Exclude cells 1, 4, 7 in the diagram).
          euyb: if (y < AuxWHD(2)) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx + FldSlv % whd(1))
          end if euyb

          NewP = (FldSlv % prhs(idx) - OffDiag) / Diag

          Err = max(Err, abs(FldSlv % p(idx) - NewP))

          FldSlv % p(idx) = NewP
          idx = idx + 1
        end do lx
      end do ly

      if (Err < Derr) return
    end do iterl

  end subroutine GSPSlv2D
!==============================================================================!

!==============================================================================!
  pure subroutine GSPSlv3D(FldSlv, MaxIt, DErr)
    ! Gauss-Seidel pressure solve.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)     , intent(in)    :: DErr   ! Desired error threshold.
    integer      , intent(in)    :: MaxIt  ! Maximum number of iterations.
    real(dp)                     :: scale, Err ! Scaling factor, iteration error.
    real(dp)                     :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                     :: NewP ! New pressure value.
    integer                      :: AuxWHD(3)
    integer                      :: iter, idx, z, y, x ! Iteration, counters, index


    ! Set scale.
    scale  = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
    !dir$ vector aligned
    AuxWHD = FldSlv % whd - 1

    ! ITERation Loop.
    !dir$ vector aligned
    iterl: do iter = 1, MaxIt
      Err = 0._dp
      idx = 0
      lz: do z = 0, AuxWHD(3)
        ! Loop over y.
        ly: do y = 0, AuxWHD(2)
          ! Loop over x.
          lx: do x = 0, AuxWHD(1)
            Diag    = 0._dp
            OffDiag = 0._dp

            ! Building the matrix as a five-point stencil. Solid borders, i.e. no fluid is outside of the simulation bounds.

            ! The pressure contribution for each cell depends on its 6 nearest othogonal neighbours.
            ! For example, cell a22 will receive contributions from cells a12, a23, a33 and a21.
            ! Boundary cells will only receive contributions from the orthogonally placed cells that are found within the simulation bounds. For example, cell a131 will only receive contributions from cells a121, a231 and a132, while cell a212 from cells a112, a222, a312, a211 and a213. This leads to the series of if statements found herein.

            ! y-axis
            ! ^   z-axis
            ! | /
            ! |__> x-axis
            !
            !                                 3 x 3 x 3 Cube
            !     Front plane                  Middle plane                 Back plane
            ! ----------------------      ----------------------      ----------------------
            ! | a131 | a231 | a331 |      | a132 | a232 | a332 |      | a133 | a233 | a333 |
            ! ----------------------      ----------------------      ----------------------
            ! | a121 | a221 | a321 |      | a122 | a222 | a322 |      | a123 | a223 | a323 |
            ! ----------------------      ----------------------      ----------------------
            ! | a111 | a211 | a311 |      | a112 | a212 | a312 |      | a113 | a213 | a313 |
            ! ----------------------      ----------------------      ----------------------
            !                                        |
            !                                        |
            !                                        | Mapped to
            !                                        |
            !                                        v
            ! [a111,a211,a311,a121,a221,a321,a131,a231,a331,a112,a212,a312,a122,a222,a322,a132,a232,a332,a113,a213,a313,a123,a223,a323,a133,a233,a333]

            ! Exclude Lower X Bound (Exclude cells a11, a12, a13 in the diagram).
            elxb: if (x > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - 1)
            end if elxb

            ! Exclude Lower X Bound (Exclude cells a13, a23, a33 in the diagram).
            euxb: if (x < AuxWHD(1)) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx + 1)
            end if euxb

            ! Exclude Upper Y Bound (Exclude cells a11, a12, 9 in the diagram).
            elyb: if (y > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1))
            end if elyb

            ! Exclude Lower Y Bound (Exclude cells 1, 4, 7 in the diagram).
            euyb: if (y < AuxWHD(2)) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx + FldSlv % whd(1))
            end if euyb

            elzb: if (z > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1) * FldSlv % whd(2))
            end if elzb

            euzb: if (z < AuxWHD(3)) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx + FldSlv % whd(1) * FldSlv % whd(2))
            end if euzb

            NewP = (FldSlv % prhs(idx) - OffDiag) / Diag

            Err = max(Err, abs(FldSlv % p(idx) - NewP))

            FldSlv % p(idx) = NewP
            idx = idx + 1
          end do lx
        end do ly
      end do lz

      if (Err < Derr) return
    end do iterl
  end subroutine GSPSlv3D
!==============================================================================!

!==============================================================================!
! Tested: 1D, 2D, 3D.
  pure function Lerp1D(FldQty, yi, x)
    ! Linear intERPolate in 1D.
    ! 1D interpolation between yi(1) < y < yi(4) for 0 < x < 1
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)       , intent(in) :: yi(2), x
    real(dp)                    :: lerp1d

    Lerp1D = yi(1)*(1._dp - x) + yi(2)*x
  end function Lerp1D
!==============================================================================!

!==============================================================================!
  elemental function Lerp2D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: VtxVal(4) ! VALues at cell VERtices. VtxVal(1) = x00, VtxVal(2) = x10, VtxVal(3) = x01, VtxVal(4) = x11
    type(Vec2D)    , intent(in) :: xi ! xi%x(1) = x, xi%x(2) = y
    real(dp)                    :: x(2)
    real(dp)                    :: Lerp2D
    integer, parameter          :: ivtc(8) = [0,1,0,1,0,0,1,1] ! Index array of the cell vertex displacements.
    integer                     :: ix(2)

    ! x01   -----------   x11
    !       |         |
    !       |         |
    !       |         |
    ! x00   -----------   x10

    ! Clamp coordinates to the boundaries.
    !dir$ vector aligned
    x = min( max(xi % x - FldQty % ost(1:2), 0._dp), real(FldQty % whd(1:2) - 1.001_dp, dp) )

    ! Prepare interpolation values.
    ix = aint( x )

    ! We make sure 0 < x < 1
    x  = x - ix

    ! Values at vertices.
    VtxVal = FldQty % ReadVal(ix(1) + ivtc(1:4), ix(2) + ivtc(5:8))

    Lerp2D = FldQty % Lerp(                        &
                [                                  &
                  ! Lower edge: x00, x10.
                  FldQty % Lerp(VtxVal(1:2),x(1)), &
                  ! Upper edge: x01, x11.
                  FldQty % Lerp(VtxVal(3:4),x(1))  &
                ],                                 &
                x(2)                               &
                          )
  end function Lerp2D
!==============================================================================!

!==============================================================================!
  elemental function Lerp3D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: VtxVal(8) ! Values at cell vertices. VtxVal(1) = x000, VtxVal(2) = x100, VtxVal(3) = x010, VtxVal(4) = x110, VtxVal(5) = x001, VtxVal(6) = x101, VtxVal(7) = x011, VtxVal(8) = x111
    type(Vec3D),     intent(in) :: xi! xi%x(1) = x, xi%x(2) = y, xi%x(3) = z
    real(dp)                    :: x(3)
    real(dp)                    :: Lerp3D
    integer, parameter          :: ivtc(24) = [0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1] ! Index array of the cell vertex displacements.
    integer                     :: ix(3)

    !    x011   ------------- x111
    !          /|          /|
    !         / |         / |
    !        /  |        /  |
    ! x010  -------------   | x110
    !       |   |       |   |
    !  x001 |   --------|--- x101
    !       |  /        |  /
    !       | /         | /
    !       |/          |/
    ! x000  ------------- x100

    ! Clamp coordinates to the boundaries.
    !dir$ vector aligned
    x = min( max(xi % x - FldQty % ost(1:3), 0._dp), real(FldQty % whd(1:3) - 1.001_dp, dp) )

    ! Prepare interpolation values.
    ix = aint( x )
    ! We make sure 0 < x < 1
    x  = x - ix

    VtxVal = FldQty % ReadVal(ix(1) + ivtc(1:8), ix(2) + ivtc(9:16), ix(3) + ivtc(17:24))

    Lerp3D = FldQty % Lerp(                            &
                [                                      &
                  ! Front face x000, x100, x010, x110
                  FldQty % Lerp(                       &
                    [                                  &
                      ! Lower edge: x000, x100.
                      FldQty % Lerp(VtxVal(1:2),x(1)), &
                      ! Upper edge: x010, x110.
                      FldQty % Lerp(VtxVal(3:4),x(1))  &
                    ]                                , &
                    x(2)                               &
                               ),                      &
                  ! Back face x001, x101, x011, x111
                  FldQty % Lerp(                       &
                    [                                  &
                      ! Lower edge: x001, x101.
                      FldQty % Lerp(VtxVal(5:6),x(1)), &
                      ! Upper edge: x011, x111.
                      FldQty % Lerp(VtxVal(7:8),x(1))  &
                    ]                                , &
                    x(2)                               &
                              )                        &
                ],                                     &
                x(3)                                   &
                          )
  end function Lerp3D
!==============================================================================!

!==============================================================================!
  pure subroutine MaxTimestep2D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, SRVel ! Eucledian norm of the velocity, Max(Vel)
    real(dp)                       :: cy, cx ! Center i, center j
    integer                        :: y, x   ! Counters

    MaxVel = 0._dp

    ly: do y = 0, FldSlv % whd(2) - 1
      ! Half way to the next y-index.
      cy = y + 0.5_dp
      lx: do x = 0, FldSlv % whd(1) - 1
        ! Half way to the next x-index.
        cx = x + 0.5_dp
        ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
        SRVel =  sum( FldSlv % u % Lerp(Vec2D([cx, cy])) * FldSlv % u % Lerp(Vec2D([cx, cy])) )
        ! Compare to Maximal Velocity.
        if (SRVel < MaxVel) cycle lx
        ! Update MaxVel
        MaxVel = SRVel
      end do lx
    end do ly

    ! Fluid should not flow through more than one cell per iteration.
    FldSlv % TimeStep(1) = FldSlv % safety * FldSlv % celsiz / sqrt(MaxVel)

    ! Clamp Time Step to be within Reasonable Bounds.
    ctsrb: if (FldSlv % TimeStep(1) < FldSlv % TimeStep(2) .or. t + FldSlv % TimeStep(1) == t) then
      ! Increase Timestep.
      it: do
        FldSlv % TimeStep(1) = FldSlv % TimeStep(1) + FldSlv % TimeStep(2)
        if (t + FldSlv % TimeStep(1) > t) return
      end do it
    else if (FldSlv % TimeStep(1) > FldSlv % TimeStep(3)) then ctsrb
      FldSlv % TimeStep(1) = FldSlv % TimeStep(3)
    end if ctsrb
  end subroutine MaxTimestep2D

  pure subroutine MaxTimestep3D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, SRVel ! Eucledian norm of the velocity, Max(Vel).
    real(dp)                       :: cz, cy, cx ! Center i, center j, center k
    integer                        :: z, y, x    ! Counters

    MaxVel = 0._dp

    lz: do z = 0, FldSlv % whd(3) - 1
      cz = z + 0.5_dp
      ly: do y = 0, FldSlv % whd(2) - 1
        ! Half way to the next y-index.
        cy = y + 0.5_dp
        lx: do x = 0, FldSlv % whd(1) - 1
          ! Half way to the next x-index.
          cx = x + 0.5_dp
          ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
          SRVel = sum( FldSlv % u % Lerp(Vec3D([cx, cy, cz])) * FldSlv % u % Lerp(Vec3D([cx, cy, cz])) )
          ! Compare to Maximal Velocity.
          cmv: if (SRVel > MaxVel) then
            ! Update MaxVel
            MaxVel = SRVel
          else cmv
            cycle lx
          end if cmv
        end do lx
      end do ly
    end do lz

    ! Fluid should not flow through more than one cell per iteration.
    FldSlv % TimeStep(1) = FldSlv % safety * FldSlv % celsiz / sqrt(MaxVel)

    ! Clamp Time Step to be within Reasonable Bounds.
    ctsrb: if (FldSlv % TimeStep(1) < FldSlv % TimeStep(2) .or. t + FldSlv % TimeStep(1) == t) then
      ! Increase Timestep.
      it: do
        FldSlv % TimeStep(1) = FldSlv % TimeStep(1) + FldSlv % TimeStep(2)
        if (t + FldSlv % TimeStep(1) > t) return
      end do it
    else if (FldSlv % TimeStep(1) > FldSlv % TimeStep(3)) then ctsrb
      FldSlv % TimeStep(1) = FldSlv % TimeStep(3)
    end if ctsrb
  end subroutine MaxTimestep3D

  pure subroutine FEulerVel2D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(FluidQty) , intent(in)  :: FldVel(2)
    class(FluidQty), intent(in)  :: FldQty
    type(Vec2D)    , intent(in)  :: xi
    real(dp)       , intent(out) :: xf(2)
    real(dp)                     :: vel(2) ! Velocities.
    real(dp)       , intent(in)  :: h
    !DIR$ ASSUME_ALIGNED xf: 8

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    !dir$ vector aligned
    xf = xi % x - vel * h
  end subroutine FEulerVel2D

  pure subroutine FEulerVel3D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(FluidQty) , intent(in)  :: FldVel(3)
    class(FluidQty), intent(in)  :: FldQty
    type(Vec3D)    , intent(in)  :: xi
    real(dp)       , intent(out) :: xf(3)
    real(dp)                     :: vel(3) ! Velocities.
    real(dp)       , intent(in)  :: h
    !DIR$ ASSUME_ALIGNED xf: 8

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    !dir$ vector aligned
    xf = xi % x - vel * h
  end subroutine FEulerVel3D

  pure subroutine Advect2D(FldQty, FldVel, h)
    implicit none
    type(FluidQty2D), intent(in)    :: FldVel
    class(FluidQty) , intent(inout) :: FldQty
    real(dp)                        :: xi(2), xf(2)
    real(dp)        , intent(in)    :: h
    integer                         :: idx, y, x

    ! Loop through y.
    idx = 0
    ly: do y = 0, FldQty % whd(2) - 1
      xi(2) = y + FldQty % ost(2)
      ! Loop through x.
      lx: do x = 0, FldQty % whd(1) - 1
        xi(1) = x + FldQty % ost(1)
        ! Integrate FldQty.
        call FldQty % FEuler(FldVel % FldQty(:), Vec2D(xi), xf, h)
        ! Interpolate component from grid and write new value.
        FldQty % new(idx) = FldQty % Lerp(Vec2D(xf))
        idx = idx + 1
      end do lx
    end do ly
  end subroutine Advect2D

  subroutine Advect3D(FldQty, FldVel, h)
    implicit none
    type(FluidQty3D), intent(in)    :: FldVel
    class(FluidQty) , intent(inout) :: FldQty
    real(dp)                        :: xi(3), xf(3)
    real(dp)        , intent(in)    :: h
    integer                         :: idx, z, y, x

    ! Loop through z.
    idx = 0
    lz: do z = 0, FldQty % whd(3) - 1
      xi(3) = z + FldQty % ost(3)
      ! Loop through y.
      ly: do y = 0, FldQty % whd(2) - 1
        xi(2) = y + FldQty % ost(2)
        ! Loop through x.
        lx: do x = 0, FldQty % whd(1) - 1
          xi(1) = x + FldQty % ost(1)
          ! Integrate FldQty.
          call FldQty % FEuler(FldVel % FldQty(:), Vec3D(xi), xf, h)
          ! Interpolate component from grid and save the new value.
          FldQty % new(idx) = FldQty % Lerp(Vec3D(xf))
          idx = idx + 1
        end do lx
      end do ly
    end do lz
  end subroutine Advect3D

  subroutine ApplyPressure2D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale, sxp ! Scale, scale x pressure.
    integer                        :: AuxWHD(2)
    integer                        :: idx, y, yp1, x, xp1

    scale = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz)

    !dir$ vector aligned
    AuxWHD = FldSlv % whd - 1

    ! Applies pressure and updates the velocity field.
    idx = 0
    !dir$ vector aligned
    ly: do y = 0, AuxWHD(2)
      yp1 = y + 1
      lx: do x = 0, AuxWHD(1)
        xp1 = x + 1
        sxp = scale * FldSlv % p(idx)
        call FldSlv % u % WriteOldVal(x, y, FldSlv % u % ReadVal(x, y) - sxp)
        call FldSlv % u(1) % WriteOldVal(xp1, y  , FldSlv % u(1) % ReadVal(xp1, y  ) + sxp)
        call FldSlv % u(2) % WriteOldVal(x  , yp1, FldSlv % u(2) % ReadVal(x  , yp1) + sxp)
        idx = idx + 1
      end do lx
    end do ly
    ! Set boundary velocities to zero (conservation of mass, nothing comes in or goes out).

    !                               at ix \in [1, width] and iy = height
    !                                              Vy = 0
    !
    !                                           ^     ^     ^
    !                                           |     |     |
    !                             (x0,y1) ------------------------- (x1,y1)
    !                                     |     |     |     |     |
    !                                     |     v     v     v     |
    !                                     |                       |
    !at ix = 1 and iy \in [1, height] <---|--->               <---|--->  at ix = width and iy \in [1, height]
    !            Vx = 0                   |                       |                   Vx = 0
    !                                 <---|--->               <---|--->
    !                                     |                       |
    !                                     |     ^     ^     ^     |
    !                                     |     |     |     |     |
    !                             (x0,y0) ------------------------- (x1,y0)
    !                                           |     |     |
    !                                           v     v     v
    !
    !                                  at ix \in [1, width] and iy = 1
    !                                              Vy = 0

    ! Vertical boundaries.
    !dir$ parallel
    !dir$ loop count min(128)
    !dir$ vector aligned
    vb: do y = 0, AuxWHD(2)
      ! Lower x boundary.
      call FldSlv % u(1) % WriteOldVal(0, y, 0._dp)
      ! Upper x boundary
      call FldSlv % u(1) % WriteOldVal(FldSlv % whd(1), y, 0._dp)
    end do vb

    ! Horizontal boundaries.
    !dir$ parallel
    !dir$ loop count min(128)
    !dir$ vector aligned
    hb: do x = 0, AuxWHD(1)
      ! Lower y boundary.
      call FldSlv % u(2) % WriteOldVal(x, 0, 0._dp)
      ! Upper y boundary
      call FldSlv % u(2) % WriteOldVal(x, FldSlv % whd(2), 0._dp)
    end do hb

  end subroutine ApplyPressure2D

  subroutine ApplyPressure3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale, sxp ! Scale, scale x pressure.
    integer                        :: AuxWHD(3)
    integer                        :: idx, z, y, x, zp1, yp1, xp1

    scale = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz)

    !dir$ vector aligned
    AuxWHD = FldSlv % whd - 1

    ! Applies pressure and updates the velocity field.
    idx = 0
    !dir$ vector aligned
    lz: do z = 0, AuxWHD(3)
      zp1 = z + 1
      ly: do y = 0, AuxWHD(2)
        yp1 = y + 1
        lx: do x = 0, AuxWHD(1)
          xp1 = x + 1
          sxp = scale * FldSlv % p(idx)
          ! Update velocity field.
          call FldSlv % u % WriteOldVal(x, y, z, FldSlv % u % ReadVal(x, y, z) - sxp)
          call FldSlv % u(1) % WriteOldVal(xp1, y  , z  , FldSlv % u(1) % ReadVal(xp1, y  , z  ) + sxp)
          call FldSlv % u(2) % WriteOldVal(x  , yp1, z  , FldSlv % u(2) % ReadVal(x  , yp1, z  ) + sxp)
          call FldSlv % u(3) % WriteOldVal(x  , y  , zp1, FldSlv % u(3) % ReadVal(x  , y  , zp1) + sxp)
          idx = idx + 1
        end do lx
      end do ly
    end do lz

    ! Set boundary velocities to zero (conservation of mass, nothing comes in or goes out).
    !        (x0,y1,z1)      -------------------------- (x1,y1,z1)
    !                       /|                       /|
    !                      / |        ^             / |
    !                     /  |        | Vy = 0     /  |
    !                    /   |        |    Vz = 0 /   |
    !                   /    |        |     ↗    /    |
    !                  /     |        |    /    /     |
    !     (x0,y1,z0)  --------------------------      | (x1,y1,z0)
    !                 |      |        v  /     |      |
    !          Vx = 0 |      |          /      |      |  Vx = 0
    !      <----------| ---------->    ↙  <----|----- ----------->
    !                 |      |                 |      |
    !                 |      |       ↗         |      |
    !                 |      |      /          |      |
    !      (x0,y0,z1) |      ------/-----------|------ (x1,y0,z1)
    !                 |     /     /            |     /
    !                 |    /     ↙   ^         |    /
    !                 |   /  Vz = 0  | Vy = 0  |   /
    !                 |  /           |         |  /
    !                 | /                      | /
    !                 |/             |         |/
    !     (x0,y0,z0)  -------------------------- (x1,y0,z0)
    !                                |
    !                                v

    !          XY Front Plane                         XY Back Plane
    ! (x0,y1,z0) ----------- (x1,y1,z0)     (x0,y1,z1) ----------- (x1,y1,z1)
    !            |         |                           |         |
    !            |         |                           |         |
    !            |         |                           |         |
    ! (x0,y0,z0) ----------- (x1,y0,z0)     (x0,y0,z1) ----------- (x1,y0,z1)
    ! Perpendicular velocity is Vz = u(3).
    ! Front Plane: At iz = 1, x \in [1, width] and y \in [1, height], Vz = 0
    ! Back Plane: At iz = depth, x \in [1, width] and y \in [1, height], Vz = 0

    !           YZ Left Plane                        YZ Right Plane
    ! (x0,y1,z1) ----------- (x0,y1,z0)     (x1,y1,z1) ----------- (x1,y1,z0)
    !            |         |                           |         |
    !            |         |                           |         |
    !            |         |                           |         |
    ! (x0,y0,z1) ----------- (x0,y0,z0)     (x1,y0,z1) ----------- (x1,y0,z0)
    ! Perpendicular velocity is Vx = u(1).
    ! Front Plane: At ix = 1, z \in [1, depth] and y \in [1, height], Vx = 0
    ! Back Plane: At ix = width, z \in [1, depth] and y \in [1, height], Vx = 0

    !          XZ Bottom Plane                        XZ Top Plane
    ! (x0,y0,z1) ----------- (x1,y0,z1)     (x0,y1,z1) ----------- (x1,y1,z1)
    !            |         |                           |         |
    !            |         |                           |         |
    !            |         |                           |         |
    ! (x0,y0,z0) ----------- (x1,y0,z0)     (x0,y1,z0) ----------- (x1,y1,z0)
    ! Perpendicular velocity is Vy = u(2).
    ! Front Plane: At iy = 1, z \in [1, depth] and x \in [1, widht], Vy = 0
    ! Back Plane: At iy = height, z \in [1, depth] and x \in [1, widht], Vy = 0

    ! YZ Plane Boundaries.
    ! Loop through z.
    !dir$ parallel
    !dir$ loop count min(4)
    !dir$ vector aligned
    lz2: do z = 0, AuxWHD(3)
      ! Loop through y.
      !dir$ ivdep
      !dir$ vector aligned
      ly3: do y = 0, AuxWHD(2)
        ! Left YZ plane (ix = 1)
        call FldSlv % u(1) % WriteOldVal(0, y, z, 0._dp)
        ! Right YZ plane (ix = whd(1))
        call FldSlv % u(1) % WriteOldVal(FldSlv % whd(1), y, z, 0._dp)
      end do ly3
    end do lz2

    ! XY Plane Boundaries
    ! Loop through y.
    !dir$ parallel
    !dir$ loop count min(8)
    !dir$ vector aligned
    ly2: do y = 0, AuxWHD(2)
      ! Loop through x.
      !dir$ ivdep
      !dir$ vector aligned
      lx2: do x = 0, AuxWHD(1)
        ! Front XY plane (iz = 1)
        call FldSlv % u(3) % WriteOldVal(x, y, 0, 0._dp)
        ! Back XY plane (iz = whd(3))
        call FldSlv % u(3) % WriteOldVal(x, y, FldSlv % whd(3), 0._dp)
      end do lx2
    end do ly2

    ! XZ Plane Boundaries.
    ! Loop through z.
    !dir$ parallel
    !dir$ loop count min(4)
    !dir$ vector aligned
    lz3: do z = 0, AuxWHD(3)
      ! Loop through x.
      !dir$ ivdep
      !dir$ vector aligned
      lx3: do x = 0, AuxWHD(1)
        ! Bottom XZ plane (iy = 1)
        call FldSlv % u(2) % WriteOldVal(x, 0, z, 0._dp)
        ! Top XZ plane (iy = whd(2))
        call FldSlv % u(2) % WriteOldVal(x, FldSlv % whd(2), z, 0._dp)
      end do lx3
    end do lz3
  end subroutine ApplyPressure3D

  pure function ReadArray2D(FldSlv, array, i, j)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    real(dp)                    :: ReadArray2D
    integer,         intent(in) :: i, j

    !dir$ vector aligned
    ReadArray2D = array(i + FldSlv % whd(1)*j)
  end function ReadArray2D

  pure function ReadArray3D(FldSlv, array, i, j, k)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    real(dp)                    :: ReadArray3D
    integer,         intent(in) :: i, j, k

    !dir$ vector aligned
    ReadArray3D = array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*FldSlv % whd(2)*k)
  end function ReadArray3D

  pure subroutine WriteArray2D(FldSlv, array, i, j, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    real(dp),        intent(in) :: value
    integer,         intent(in) :: i, j

    !dir$ vector aligned
    array(i + FldSlv % whd(1)*j) = value
  end subroutine WriteArray2D

  pure subroutine WriteArray3D(FldSlv, array, i, j, k, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    real(dp),        intent(in) :: value
    integer,         intent(in) :: i, j, k

    !dir$ vector aligned
    array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*FldSlv % whd(2)*k) = value
  end subroutine WriteArray3D

  elemental function ReadVal2D(FldQty, i, j)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal2D
    integer,         intent(in) :: i, j

    !dir$ vector aligned
    ReadVal2D = FldQty % old(i + FldQty % whd(1)*j)
  end function ReadVal2D

  elemental function ReadVal3D(FldQty, i, j, k)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal3D
    integer,         intent(in) :: i, j, k

    !dir$ vector aligned
    ReadVal3D = FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1)*FldQty % whd(2)*k)
  end function ReadVal3D

  elemental subroutine WriteNewVal2D(FldQty, i, j, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j

    !dir$ vector aligned
    FldQty % new(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteNewVal2D

  elemental subroutine WriteNewVal3D(FldQty, i, j, k, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j, k

    !dir$ vector aligned
    FldQty % new(i + FldQty % whd(1)*j + FldQty % whd(1) * FldQty % whd(2) * k) = NewVal
  end subroutine WriteNewVal3D

  elemental subroutine WriteOldVal2D(FldQty, i, j, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j

    !dir$ vector aligned
    FldQty % old(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteOldVal2D

  elemental subroutine WriteOldVal3D(FldQty, i, j, k, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j, k

    !dir$ vector aligned
    FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1)*FldQty % whd(2)*k) = NewVal
  end subroutine WriteOldVal3D

  elemental subroutine UpdateVals(FldQty)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    !dir$ loop count min(256)
    !dir$ vector aligned
    FldQty % old = FldQty % new
    !dir$ loop count min(512)
    !dir$ vector aligned
    FldQty % new = 0.
  end subroutine UpdateVals

  subroutine WriteToFile(FldSlv, unit)
    class(FluidSlv), intent(in)   :: FldSlv
    integer, intent(in), optional :: unit
    integer, allocatable, save    :: ndim
    integer :: i

    if(.not. allocated(ndim)) allocate(ndim, source = size(FldSlv % whd))

    ! Check if Unit is Present.
    cup: if (present(unit)) then
      ! Select Unit.
      su: select case (unit)
        case (1) su
          ! Print density.
          open(unit=unit, file="density.dat", status="replace", access="stream")
          !dir$ vector aligned
          write(unit) FldSlv % rho % old
          close(unit=unit)
        case (2) su
          ! Print velocities.
          open(unit=unit, file="velocity.dat", status="replace", access="stream")
          ! Loop through the Number of DIMensions.
          lndim: do i = 1, ndim
            !dir$ vector aligned
            write(unit) FldSlv % u(i) % old
          end do lndim
          close(unit=unit)
        case (3) su
          ! Print pressures.
          open(unit=unit, file="pressure.dat", status="replace", access="stream")
          !dir$ vector aligned
          write(unit) FldSlv % p
          close(unit=unit)
      end select su
    else cup
      ! If unit is not present. Print everything.
      open(unit=1, file="density.dat", status="replace", access="stream")
      !dir$ vector aligned
      write(1) FldSlv % rho % old
      close(unit=1)

      open(unit=2, file="velocity.dat", status="replace", access="stream")
      ! Loop through the Number of DIMensions.
      lndim2: do i = 1, ndim
        !dir$ vector aligned
        write(2) FldSlv % u(i) % old
      end do lndim2
      close(unit=2)

      open(unit=3, file="pressure.dat", status="replace", access="stream")
      write(3) FldSlv % p
      !dir$ vector aligned
      close(unit=3)
    end if cup
  end subroutine WriteToFile

end module fluid
