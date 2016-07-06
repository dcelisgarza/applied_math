! Generalised Fortran version of original c++ code by TunnaBrain (Benedikt Bitterli) https://github.com/tunabrain/incremental-fluids.

module fluid
  use nrtype

  type FluidQty
    character(len=:), allocatable :: name   ! Variable name.
    real(dp)        , allocatable :: old(:) ! Old value.
    real(dp)        , allocatable :: new(:) ! New value.
    integer         , allocatable :: whd(:) ! Width, height, depth of the simulation domain.
    real(dp)        , allocatable :: ost(:) ! Offset.
    real(dp)                      :: celsiz ! Cell size (area or volume).
  contains
    ! Initialise type.
    procedure :: Init => InitFluidQty
    ! Add fluid flow.
    procedure :: AddFlow2D
    procedure :: AddFlow3D
    generic   :: AddFlow => AddFlow2D, AddFlow3D
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
    integer        , allocatable :: whd(:)   ! Width, height, depth
    real(dp)                     :: celsiz   ! Cell size.
    real(dp)                     :: srho     ! Scalar (fluid) density.
    real(dp)       , allocatable :: prhs(:)  ! Right hand side of pressure equation.
    real(dp)       , allocatable :: p(:)     ! Pressure solution.
    real(dp)                     :: TimeStep, MaxTimeStep, MinTimeStep ! Timesteps
    real(dp)                     :: safety   ! Safety factor for timestep calculation
  contains
    ! Initialise type.
    procedure :: Init => InitFluid
    ! Add flow.
    procedure :: InitFlow2D
    procedure :: InitFlow3D
    generic   :: AddFlow => InitFlow2D, InitFlow3D
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
    procedure :: ApplyPressure2D, ApplyPressure3D
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

  type Vec2D
    real(dp) :: x(2)
  end type Vec2D

  type Vec2D2
    real(dp) :: x(4)
  end type Vec2D2

  type FluidQty2D
    type(FluidQty) :: FldQty(2)
  end type FluidQty2D

  type Vec3D
    real(dp) :: x(3)
  end type Vec3D

  type Vec3D2
    real(dp) :: x(6)
  end type Vec3D2

  type FluidQty3D
    type(FluidQty) :: FldQty(3)
  end type FluidQty3D

contains

!==============================================================================!
! Tested.
  subroutine InitFluid(FldSlv, names, whd, rho, timestep, scale, safety)
    implicit none
    class(FluidSlv) , intent(out)  :: FldSlv ! Fluid
    character(len=*), dimension(:), intent(in) :: names
    integer         , intent(in)   :: whd(:) ! Width (i), height (j), depth (k)
    real(dp)        , intent(in)   :: rho
    real(dp)        , intent(in)   :: timestep
    real(dp), optional, intent(in) :: scale
    real(dp), optional, intent(in) :: safety
    real(dp), dimension(size(whd)) :: ost    ! Offset
    real(dp), dimension(size(whd)) :: TMPost
    integer                        :: ndim   ! Number of dimensions
    integer                        :: ncells ! Number of cells
    integer, dimension(size(whd))  :: TMncells ! Number of cells for each different quantity
    integer                        :: i      ! Counter

    ! Number of dimensions.
    ndim = size(whd)
    ncells = product(whd) - 1

    ! Set the fluid density.
    FldSlv % srho = rho

    ! Set initial time step.
    FldSlv % TimeStep = timestep
    if(present(safety) .eqv. .true.) then
      FldSlv % MinTimeStep = timestep * safety * 1d-3
      FldSlv % MaxTimeStep = timestep * safety * 1d+3
    else
      FldSlv % MinTimeStep = timestep * 0.8_dp * 1d-3
      FldSlv % MaxTimeStep = timestep * 0.8_dp * 1d+3
    end if


    ! Cell size.
    if(present(scale) .eqv. .true.) then
      FldSlv % celsiz = scale/minval(whd)
    else
      FldSlv % celsiz = 1./minval(whd)
    endif

    ! Safety factor.
    if(present(safety) .eqv. .true.) then
      FldSlv % safety = safety
    else
      FldSlv % safety = 0.8_dp
    endif

    ! Initialise density.
    ! Centered in the grid.
    ost    = 0.5
    call FldSlv % rho % init(trim(names(1)(1:len(names(1)))), whd, ost, FldSlv % celsiz)

    ! Initialise velocities.
    ! Staggered on the grid.
    allocate(FldSlv % u(ndim))
    iv: do i = 1, ndim
      TMncells    = whd
      TMncells(i) = whd(i) + 1
      TMPost    = ost
      TMPost(i) = ost(i) - 0.5
      call FldSlv % u(i) % init(trim(names(i+1)(1:len(names(i+1)))), TMncells, TMPost, FldSlv % celsiz)
    end do iv

    ! Allocate whd (simulation bounds).
    allocate (FldSlv % whd(ndim))
    FldSlv % whd = whd

    ! Allocate the right hand side of pressure solve.
    allocate(FldSlv % prhs(0:ncells))
    !dir$ loop count min(256)
    FldSlv % prhs = 0._dp

    ! Allocate the pressure solution and give initial guess.
    allocate(FldSlv % p(0:ncells))
    !dir$ loop count min(256)
    FldSlv % p = 0._dp
  end subroutine InitFluid
!==============================================================================!

!==============================================================================!
! Tested
  subroutine InitFluidQty(FldQty, name, whd, ost, celsiz)
    implicit none
    character(len=*), intent(in)  :: name   ! Variable name.
    integer,          intent(in)  :: whd(:) ! Width (i), height (j), depth (k)
    real(dp),         intent(in)  :: ost(:) ! Offset
    real(dp),         intent(in)  :: celsiz ! Cell size (area or volume)
    class(FluidQty),  intent(out) :: FldQty ! Fluid Quantity
    integer                       :: ncells ! # of cells = prod(whd)

    ! Check that size(whd) == size(ost).
    check_size: if ( size(whd) /= size(ost) ) then
      write(*,*) " Error: fluid.f08: InitFluidQty: check_size: Variable name = ", name, "."
      write(*,*) " Size(whd) == size(ost); size(whd) = ", size(whd), " size(ost) = ", size(ost)
      stop
    end if check_size

    ! Allocate name.
    allocate(character(len=len(name)):: FldQty % name)
    FldQty % name = name

    ! Width, height, depth
    allocate(FldQty % whd(size(whd)))
    FldQty % whd = whd

    ! Offset
    allocate(FldQty % ost(size(ost)))
    FldQty % ost = ost

    ! Allocate values.
    ncells = product(whd) - 1
    allocate(FldQty % old(0:ncells), FldQty % new(0:ncells))
    !dir$ loop count min(512)
    FldQty % old = 0.
    !dir$ loop count min(512)
    FldQty % new = 0.

    ! Cell size
    FldQty % celsiz = celsiz
  end subroutine InitFluidQty
!==============================================================================!

!==============================================================================!
  ! Initialise flow.
  ! Adds the respective quantities inside the box bounded by the coordinates xi%x.
  ! Tested: 2D
  pure subroutine InitFlow2D(FldSlv, xi, values)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    type(Vec2D2)   , intent(in)    :: xi        ! Coordinates for the vertices of the 2D box.
    real(dp)       , intent(in)    :: values(3) ! value(1) := density, value(2) = Vx, value(3) = Vy.

    call FldSlv % rho    % AddFlow(xi, values(1)  )
    call FldSlv % u(1:2) % AddFlow(xi, values(2:3))
  end subroutine InitFlow2D

  pure subroutine InitFlow3D(FldSlv, xi, values)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    type(Vec3D2)   , intent(in)    :: xi        ! Coordinates for the vertices of the 3D box.
    real(dp)       , intent(in)    :: values(4) ! value(1) := density, value(2) = Vx, value(3) = Vy, value(4) = Vz

    call FldSlv % rho    % AddFlow(xi, values(1)  )
    call FldSlv % u(1:3) % AddFlow(xi, values(2:4))
  end subroutine InitFlow3D
!==============================================================================!

!==============================================================================!
  ! Add relevant quantities to the fluid simulation.
  ! Tested: 2D.
  elemental subroutine AddFlow2D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec2D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1
    real(dp)       , intent(in)    :: value
    integer                        :: i, j, ixi(4), idx

    ! (x0,y1) ----------- (x1,y1)
    !         |         |
    !         |         |
    !         |         |
    ! (x0,y0) ----------- (x1,y0)

    ! Set up the simulation bounds.
    bounds: do concurrent (i = 1: 2)
      j = 2*i
      ! Using fortran's standard indexing we need to add 1 to the coordinates if we want to properly account for zero.
      ixi(j-1 : j) = aint( (xi % x(j-1 : j) ) / FldQty % celsiz - FldQty % ost(i) )
    end do bounds

    ! Loop through y.
    ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
    ly: do concurrent (j = max(ixi(3), 0): min( ixi(4), FldQty % whd(2)) - 1)
      ! Loop through x.
      ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
      lx: do concurrent (i = max(ixi(1), 0): min( ixi(2), FldQty % whd(1)) - 1)
        ! Write value to old value.
        idx = i + FldQty % whd(1)*j
        if ( abs(FldQty % old(idx)) < abs(value) ) FldQty % old(idx) = value
      end do lx
    end do ly
  end subroutine AddFlow2D
!==============================================================================!

!==============================================================================!
  elemental subroutine AddFlow3D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec3D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1, xi%x(5) = z0, xi%x(6) = z1
    real(dp)       , intent(in)    :: value
    integer                        :: i, j, k, ixi(6), idx

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

    bounds: do concurrent (i = 1: 3)
      j = 2*i
      ixi(j-1 : j) = aint( (xi % x(j-1 : j) ) / FldQty % celsiz - FldQty % ost(i) )
    end do bounds

    ! Loop through z.
    ! Start from the maximum value between the desired z0 and 1 (for the cell index) up to the minimum value between the desired z1 and the simulation depth.
    lz: do concurrent (k = max(ixi(5), 0): min(ixi(6), FldQty % whd(3)) - 1)
      ! Loop through y.
      ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
      ly: do concurrent (j = max(ixi(3), 0): min(ixi(4), FldQty % whd(2)) - 1)
        ! Loop through x.
        ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
        lx: do concurrent (i = max(ixi(1), 0): min(ixi(2), FldQty % whd(1)) - 1)
          ! Write value to old value.
          idx = i + FldQty % whd(1)*j + FldQty % whd(2)*k
          if ( abs(FldQty % old(idx)) < abs(value) ) FldQty % old(idx) = value
        end do lx
      end do ly
    end do lz
  end subroutine AddFlow3D
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
    integer                        :: i, j, idx ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz

    idx = 0
    ly: do j = 0, FldSlv % whd(2) - 1
      lx: do i = 0, FldSlv % whd(1) - 1
        FldSlv % prhs(idx) = -scale * (                                     &
                                         FldSlv % u(1) % ReadVal(i+1, j  )  &
                                       + FldSlv % u(2) % ReadVal(i  , j+1)  &
                                       - sum(FldSlv % u % ReadVal(i , j  )) &
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
    integer                        :: i, j, k, idx ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz
    idx = 0
    lz: do k = 0, FldSlv % whd(3) - 1
      ly: do j = 0, FldSlv % whd(2) - 1
        lx: do i = 0, FldSlv % whd(1) - 1
          FldSlv % prhs(idx) = -scale * (              &
                                           FldSlv % u(1) % ReadVal(i+1, j  , k  )  &
                                         + FldSlv % u(2) % ReadVal(i  , j+1, k  )  &
                                         + FldSlv % u(3) % ReadVal(i  , j  , k+1)  &
                                         - sum(FldSlv % u % ReadVal(i , j  , k  )) &
                                        )
          idx = idx + 1
        end do lx
      end do ly
    end do lz
  end subroutine BuildRhs3D
  !==============================================================================!

  !==============================================================================!
  ! Tested: 2D.
  subroutine GSPSlv2D(FldSlv, MaxIt, DErr)
    ! Gauss-Seidel pressure solve.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    integer        , intent(in)    :: MaxIt  ! Maximum number of iterations.
    real(dp)       , intent(in)    :: DErr   ! Desired error threshold.
    real(dp)                       :: scale, Err ! Scaling factor, iteration error.
    integer                        :: iter, i, j, idx ! Iteration, counters, index
    real(dp)                       :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                       :: NewP ! New pressure value.
    integer                        :: AuxWHD(2)

    ! Set scale.
    scale  = FldSlv % TimeStep / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
    AuxWHD = FldSlv % whd - 1

    ! ITERation Loop.
    iterl: do iter = 1, MaxIt
      Err = 0._dp
      ! Loop over y.
      idx = 0
      ly: do j = 0, AuxWHD(2)
        ! Loop over x.
        lx: do i = 0, AuxWHD(1)
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
          elxb: if (i > 0) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx - 1)
          end if elxb

          ! Exclude Lower X Bound (Exclude cells a13, a23, a33 in the diagram).
          euxb: if (i < AuxWHD(1)) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx + 1)
          end if euxb

          ! Exclude Upper Y Bound (Exclude cells a11, a12, 9 in the diagram).
          elyb: if (j > 0) then
            Diag    = Diag + scale
            OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1))
          end if elyb

          ! Exclude Lower Y Bound (Exclude cells 1, 4, 7 in the diagram).
          euyb: if (j < AuxWHD(2)) then
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
  subroutine GSPSlv3D(FldSlv, MaxIt, DErr)
    ! Gauss-Seidel pressure solve.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    integer      , intent(in)    :: MaxIt  ! Maximum number of iterations.
    real(dp)     , intent(in)    :: DErr   ! Desired error threshold.
    real(dp)                     :: scale, Err ! Scaling factor, iteration error.
    integer                      :: iter, i, j, k, idx ! Iteration, counters, index
    real(dp)                     :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                     :: NewP ! New pressure value.
    integer                      :: AuxWHD(3)

    ! Set scale.
    scale  = FldSlv % TimeStep / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
    AuxWHD = FldSlv % whd - 1

    ! ITERation Loop.
    iterl: do iter = 1, MaxIt
      Err = 0._dp
      idx = 0
      lz: do k = 0, AuxWHD(3)
        ! Loop over y.
        ly: do j = 0, AuxWHD(2)
          ! Loop over x.
          lx: do i = 0, AuxWHD(1)
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
            elxb: if (i > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - 1)
            end if elxb

            ! Exclude Lower X Bound (Exclude cells a13, a23, a33 in the diagram).
            euxb: if (i < AuxWHD(1)) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx + 1)
            end if euxb

            ! Exclude Upper Y Bound (Exclude cells a11, a12, 9 in the diagram).
            elyb: if (j > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1))
            end if elyb

            ! Exclude Lower Y Bound (Exclude cells 1, 4, 7 in the diagram).
            euyb: if (j < AuxWHD(2)) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx + FldSlv % whd(1))
            end if euyb

            elzb: if (k > 0) then
              Diag    = Diag + scale
              OffDiag = OffDiag - scale * FldSlv % p(idx - FldSlv % whd(1) * FldSlv % whd(2))
            end if elzb

            euzb: if (k < AuxWHD(3)) then
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
    type(Vec2D)    , intent(in) :: xi ! xi%x(1) = x, xi%x(2) = y
    real(dp)                    :: Lerp2D
    real(dp)                    :: x(2)
    integer                     :: ix(2)
    real(dp)                    :: VtxVal(4) ! VALues at cell VERtices. VtxVal(1) = x00, VtxVal(2) = x10, VtxVal(3) = x01, VtxVal(4) = x11
    integer, parameter          :: ivtc(8) = [0,1,0,1,0,0,1,1] ! Index array of the cell vertex displacements.

    ! x01   -----------   x11
    !       |         |
    !       |         |
    !       |         |
    ! x00   -----------   x10

    ! Clamp coordinates to the boundaries.
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
    type(Vec3D),     intent(in) :: xi! xi%x(1) = x, xi%x(2) = y, xi%x(3) = z
    real(dp)                    :: Lerp3D
    real(dp)                    :: x(3)
    integer                     :: ix(3)
    real(dp)                    :: VtxVal(8) ! Values at cell vertices. VtxVal(1) = x000, VtxVal(2) = x100, VtxVal(3) = x010, VtxVal(4) = x110, VtxVal(5) = x001, VtxVal(6) = x101, VtxVal(7) = x011, VtxVal(8) = x111
    integer, parameter          :: ivtc(24) = [0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1] ! Index array of the cell vertex displacements.

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
  subroutine MaxTimestep2D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, RMSVel ! Max and Root-mean-squared (RMS) Velocity * sqrt(n) at the centre of the grid.
    integer                        :: i, j   ! Counters
    real(dp)                       :: ci, cj ! Center i, center j

    MaxVel = 0.

    ly: do j = 0, FldSlv % whd(2) - 1
      ! Half way to the next y-index.
      cj = j + 0.5
      lx: do i = 0, FldSlv % whd(1) - 1
        ! Half way to the next x-index.
        ci = i + 0.5
        ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
        RMSVel =  sum( FldSlv % u % Lerp(Vec2D([ci, cj])) * FldSlv % u % Lerp(Vec2D([ci, cj])) )
        ! Compare to Maximal Velocity.
        if (RMSVel < MaxVel) cycle lx
        ! Update MaxVel
        MaxVel = RMSVel
      end do lx
    end do ly

    ! Fluid should not flow through more than one cell per iteration.
    FldSlv % TimeStep = FldSlv % safety * FldSlv % celsiz / sqrt(MaxVel)

    ! Clamp Time Step to be within Reasonable Bounds.
    ctsrb: if (FldSlv % TimeStep < FldSlv % MinTimeStep .or. t + FldSlv % TimeStep == t) then
      ! Increase Timestep.
      it: do
        FldSlv % TimeStep = FldSlv % TimeStep + FldSlv % MinTimeStep
        if (t + FldSlv % TimeStep > t) return
      end do it
    else if (FldSlv % TimeStep > FldSlv % MaxTimeStep) then ctsrb
      FldSlv % Timestep = (FldSlv % MinTimeStep + FldSlv % MaxTimeStep)/2.
    end if ctsrb
  end subroutine MaxTimestep2D

  subroutine MaxTimestep3D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, RMSVel ! Max and Root-mean-squared (RMS) Velocity * sqrt(n) at the centre of the grid.
    integer                        :: i, j, k    ! Counters
    real(dp)                       :: ci, cj, ck ! Center i, center j, center k

    MaxVel = 0.

    lz: do k = 0, FldSlv % whd(3) - 1
      ck = k + 0.5
      ly: do j = 0, FldSlv % whd(2) - 1
        ! Half way to the next y-index.
        cj = j + 0.5
        lx: do i = 0, FldSlv % whd(1) - 1
          ! Half way to the next x-index.
          ci = i + 0.5
          ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
          RMSVel = sum( FldSlv % u % Lerp(Vec3D([ci, cj, ck])) * FldSlv % u % Lerp(Vec3D([ci, cj, ck])) )
          ! Compare to Maximal Velocity.
          cmv: if (RMSVel > MaxVel) then
            ! Update MaxVel
            MaxVel = RMSVel
          else cmv
            cycle lx
          end if cmv
        end do lx
      end do ly
    end do lz

    ! Fluid should not flow through more than one cell per iteration.
    FldSlv % TimeStep = FldSlv % safety * FldSlv % celsiz / sqrt(MaxVel)

    ! Clamp Time Step to be within Reasonable Bounds.
    ctsrb: if (FldSlv % TimeStep < FldSlv % MinTimeStep .or. t + FldSlv % TimeStep == t) then
      ! Increase Timestep.
      it: do
        FldSlv % TimeStep = FldSlv % TimeStep + FldSlv % MinTimeStep
        if (t + FldSlv % TimeStep > t) return
      end do it
    else if (FldSlv % TimeStep > FldSlv % MaxTimeStep) then ctsrb
      FldSlv % Timestep = (FldSlv % MinTimeStep + FldSlv % MaxTimeStep)/2.
    end if ctsrb
  end subroutine MaxTimestep3D

  subroutine FEulerVel2D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(Vec2D)    , intent(in)  :: xi
    real(dp)       , intent(in)  :: h
    class(FluidQty), intent(in)  :: FldQty
    type(FluidQty) , intent(in)  :: FldVel(2)
    real(dp)       , intent(out) :: xf(2)
    real(dp)                     :: vel(2) ! Velocities.

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    xf = xi % x - vel * h
  end subroutine FEulerVel2D

  subroutine FEulerVel3D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(Vec3D)    , intent(in)  :: xi
    real(dp)       , intent(in)  :: h
    class(FluidQty), intent(in)  :: FldQty
    type(FluidQty) , intent(in)  :: FldVel(3)
    real(dp)       , intent(out) :: xf(3)
    real(dp)                     :: vel(3) ! Velocities.

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    xf = xi % x - vel * h
  end subroutine FEulerVel3D

  subroutine Advect2D(FldQty, FldVel, h)
    implicit none
    class(FluidQty) , intent(inout) :: FldQty
    type(FluidQty2D), intent(in)    :: FldVel
    real(dp)        , intent(in)    :: h
    real(dp)                        :: x(2)
    integer                         :: i, j, idx

    ! Loop through y.
    idx = 0
    ly: do j = 0, FldQty % whd(2) - 1
      ! Loop through x.
      lx: do i = 0, FldQty % whd(1) - 1
        x(1) = i + FldQty % ost(1)
        x(2) = j + FldQty % ost(2)
        ! Integrate FldQty.
        call FldQty % FEuler(FldVel % FldQty(:), Vec2D(x), x, h)
        ! Interpolate component from grid and write new value.
        FldQty % new(idx) = FldQty % Lerp(Vec2D(x))
        print*, FldQty % new(idx)
        idx = idx + 1
      end do lx
    end do ly
  end subroutine Advect2D

  subroutine Advect3D(FldQty, FldVel, h)
    implicit none
    class(FluidQty) , intent(inout) :: FldQty
    type(FluidQty3D), intent(inout) :: FldVel
    real(dp)        , intent(in)    :: h
    real(dp)                        :: x(3)
    integer                         :: i, j, k

    ! Loop through z.
    lz: do k = 0, FldQty % whd(3) - 1
      ! Loop through y.
      ly: do j = 0, FldQty % whd(2) - 1
        ! Loop through x.
        lx: do i = 0, FldQty % whd(1) - 1
          x(1) = i + FldQty % ost(1)
          x(2) = j + FldQty % ost(2)
          x(3) = k + FldQty % ost(3)
          ! Integrate FldQty.
          call FldQty % FEuler(FldVel % FldQty(:), Vec3D(x), x, h)
          ! Interpolate component from grid and save the new value.
          call FldQty % WriteNewVal(i, j, k, FldQty % Lerp(Vec3D(x)))
        end do lx
      end do ly
    end do lz
  end subroutine Advect3D

  subroutine ApplyPressure2D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale, sxp ! Scale, scale x pressure.
    integer                        :: i, j, ip1, jp1, idx, idxxp1, idxyp1

    scale = FldSlv % TimeStep / (FldSlv % srho * FldSlv % celsiz)
    !print*, scale

    ! Applies pressure and updates the velocity field.
    idx = 0
    ly: do j = 0, FldSlv % whd(2) - 1
      jp1 = j + 1
      idxxp1 = FldSlv % whd(1) * j
      idxyp1 = FldSlv % whd(1) * jp1
      lx: do i = 0, FldSlv % whd(1) - 1
        ip1 = i + 1
        idxxp1 = ip1 + idxxp1
        idxyp1 = i   + idxyp1
        sxp = scale * FldSlv % p(idx)
        call FldSlv % u % WriteOldVal(i, j, FldSlv % u % ReadVal(i, j) - sxp)
        call FldSlv % u(1) % WriteOldVal(ip1, j  , FldSlv % u(1) % ReadVal(ip1, j  ) + sxp)
        call FldSlv % u(2) % WriteOldVal(i  , jp1, FldSlv % u(2) % ReadVal(i  , jp1) + sxp)
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
    vb: do j = 0, FldSlv % whd(2) - 1
      ! Lower x boundary.
      call FldSlv % u(1) % WriteOldVal(0, j, 0._dp)
      ! Upper x boundary
      call FldSlv % u(1) % WriteOldVal(FldSlv % whd(1), j, 0._dp)
    end do vb

    ! Horizontal boundaries.
    !dir$ parallel
    !dir$ loop count min(128)
    hb: do i = 0, FldSlv % whd(1) - 1
      ! Lower y boundary.
      call FldSlv % u(2) % WriteOldVal(i, 0, 0._dp)
      ! Upper y boundary
      call FldSlv % u(2) % WriteOldVal(i, FldSlv % whd(2), 0._dp)
    end do hb

  end subroutine ApplyPressure2D

  subroutine ApplyPressure3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale, sxp ! Scale, scale x pressure.
    integer                        :: i, j, k, ip1, jp1, kp1, idx

    scale = FldSlv % TimeStep / (FldSlv % srho * FldSlv % celsiz)

    ! Applies pressure and updates the velocity field.
    idx = 0
    lz: do k = 0, FldSlv % whd(3) - 1
      kp1 = k + 1
      ly: do j = 0, FldSlv % whd(2) - 1
        jp1 = j + 1
        lx: do i = 0, FldSlv % whd(1) - 1
          ip1 = i + 1
          sxp = scale * FldSlv % p(idx)
          ! Update velocity field.init
          call FldSlv % u % WriteOldVal(i, j, k, FldSlv % u % ReadVal(i, j, k) - sxp)
          call FldSlv % u(1) % WriteOldVal(ip1, j  , k  , FldSlv % u(1) % ReadVal(ip1, j  , k  ) + sxp)
          call FldSlv % u(2) % WriteOldVal(i  , jp1, k  , FldSlv % u(2) % ReadVal(i  , jp1, k  ) + sxp)
          call FldSlv % u(3) % WriteOldVal(i  , j  , kp1, FldSlv % u(3) % ReadVal(i  , j  , kp1) + sxp)
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

    ! XY Plane Boundaries
    ! Loop through y
    !dir$ parallel
    !dir$ loop count min(64)
    ly2: do j = 0, FldSlv % whd(2) - 1
      ! Loop through x.
      !dir$ parallel
      !dir$ loop count min(64)
      lx2: do i = 0, FldSlv % whd(1) - 1
        ! Front XY plane (iz = 1)
        call FldSlv % u(3) % WriteOldVal(i, j, 1, 0._dp)
        ! Back XY plane (iz = whd(3))
        call FldSlv % u(3) % WriteOldVal(i, j, FldSlv % whd(3), 0._dp)
      end do lx2
    end do ly2

    ! YZ Plane Boundaries.
    ! Loop through z.
    !dir$ parallel
    !dir$ loop count min(64)
    lz2: do k = 0, FldSlv % whd(3) - 1
      ! Loop through y.
      !dir$ parallel
      !dir$ loop count min(64)
      ly3: do j = 0, FldSlv % whd(2) - 1
        ! Left YZ plane (ix = 1)
        call FldSlv % u(1) % WriteOldVal(1, j, k, 0._dp)
        ! Right YZ plane (ix = whd(1))
        call FldSlv % u(1) % WriteOldVal(FldSlv % whd(1), j, k, 0._dp)
      end do ly3
    end do lz2

    ! XZ Plane Boundaries.
    ! Loop through z.
    !dir$ parallel
    !dir$ loop count min(4)
    lz3: do k = 0, FldSlv % whd(3) - 1
      ! Loop through x.
      !dir$ ivdep
      lx3: do i = 0, FldSlv % whd(1) - 1
        ! Bottom XZ plane (iy = 1)
        call FldSlv % u(2) % WriteOldVal(i, 1, k, 0._dp)
        ! Top XZ plane (iy = whd(2))
        call FldSlv % u(2) % WriteOldVal(i, FldSlv % whd(2), k, 0._dp)
      end do lx3
    end do lz3
  end subroutine ApplyPressure3D

  pure function ReadArray2D(FldSlv, array, i, j)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    integer,         intent(in) :: i, j
    real(dp)                    :: ReadArray2D

    ReadArray2D = array(i + FldSlv % whd(1)*j)
  end function ReadArray2D

  pure function ReadArray3D(FldSlv, array, i, j, k)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    integer,         intent(in) :: i, j, k
    real(dp)                    :: ReadArray3D

    ReadArray3D = array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*k)
  end function ReadArray3D

  pure subroutine WriteArray2D(FldSlv, array, i, j, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    integer,         intent(in) :: i, j
    real(dp),        intent(in) :: value

    array(i + FldSlv % whd(1)*j) = value
  end subroutine WriteArray2D

  pure subroutine WriteArray3D(FldSlv, array, i, j, k, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    integer,         intent(in) :: i, j, k
    real(dp),        intent(in) :: value

    array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*k) = value
  end subroutine WriteArray3D

  elemental function ReadVal2D(FldQty, i, j)
    implicit none
    integer,         intent(in) :: i, j
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal2D

    ReadVal2D = FldQty % old(i + FldQty % whd(1)*j)
  end function ReadVal2D

  elemental function ReadVal3D(FldQty, i, j, k)
    implicit none
    integer,         intent(in) :: i, j, k
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal3D

    ReadVal3D = FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1) * FldQty % whd(2)*k)
  end function ReadVal3D

  elemental subroutine WriteNewVal2D(FldQty, i, j, NewVal)
    implicit none
    integer,         intent(in)    :: i, j
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty

    FldQty % new(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteNewVal2D

  elemental subroutine WriteNewVal3D(FldQty, i, j, k, NewVal)
    implicit none
    integer,         intent(in)    :: i, j, k
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty

    FldQty % new(i + FldQty % whd(1)*j + FldQty % whd(1) * FldQty % whd(2)*k) = NewVal
  end subroutine WriteNewVal3D

  elemental subroutine WriteOldVal2D(FldQty, i, j, NewVal)
    implicit none
    integer,         intent(in)    :: i, j
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty

    FldQty % old(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteOldVal2D

  elemental subroutine WriteOldVal3D(FldQty, i, j, k, NewVal)
    implicit none
    integer,         intent(in)    :: i, j, k
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty

    FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1) * FldQty % whd(2)*k) = NewVal
  end subroutine WriteOldVal3D

  elemental subroutine UpdateVals(FldQty)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    !dir$ loop count min(256)
    FldQty % old = FldQty % new
    !dir$ loop count min(512)
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
          write(unit) FldSlv % rho % old
          close(unit=unit)
        case (2) su
          ! Print velocities.
          open(unit=unit, file="velocity.dat", status="replace", access="stream")
          ! Loop through the Number of DIMensions.
          lndim: do i = 1, ndim
            write(unit) FldSlv % u(i) % old
          end do lndim
          close(unit=unit)
        case (3) su
          ! Print pressures.
          open(unit=unit, file="pressure.dat", status="replace", access="stream")
          write(unit) FldSlv % p
          close(unit=unit)
      end select su
    else cup
      ! If unit is not present. Print everything.
      open(unit=1, file="density.dat", status="replace", access="stream")
      write(1) FldSlv % rho % old
      close(unit=1)

      open(unit=2, file="velocity.dat", status="replace", access="stream")
      ! Loop through the Number of DIMensions.
      lndim2: do i = 1, ndim
        write(2) FldSlv % u(i) % old
      end do lndim2
      close(unit=2)

      open(unit=3, file="pressure.dat", status="replace", access="stream")
      write(3) FldSlv % p
      close(unit=3)
    end if cup
  end subroutine WriteToFile

end module fluid
