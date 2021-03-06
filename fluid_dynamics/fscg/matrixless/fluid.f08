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

    ! Number of dimensions.
    ndim   = size(whd)
    ncells = product(whd) - 1
    ost    = 0.5_dp

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
      TMncells    = whd
      TMncells(i) = whd(i) + 1
      TMPost    = ost
      TMPost(i) = ost(i) - 0.5_dp
      call FldSlv % u(i) % init(trim(names(i+1)(1:len(names(i+1)))), TMncells, TMPost, FldSlv % celsiz)
    end do iv

    ! Initialise density.
    ! Centered in the grid.
    call FldSlv % rho % init(trim(names(1)(1:len(names(1)))), whd, ost, FldSlv % celsiz)

    ! Allocate the right hand side of pressure solve.
    allocate(FldSlv % prhs(0:ncells))
    FldSlv % prhs = 0._dp

    ! Allocate the pressure solution and give initial guess.
    allocate(FldSlv % p(0:ncells))
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
    FldSlv % srho = rho

    ! Safety factor.
    if(present(safety) .eqv. .true.) then
      FldSlv % safety = safety
    else
      FldSlv % safety = 0.8_dp
    endif

    ! Allocate whd (simulation bounds).
    allocate (FldSlv % whd(ndim))
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

    ! Allocate values.
    ncells = product(whd) - 1
    allocate(FldQty % old(0:ncells), FldQty % new(0:ncells))
    FldQty % old = 0._dp
    FldQty % new = 0._dp

    ! Offset
    allocate(FldQty % ost(size(ost)))
    FldQty % ost = ost

    ! Cell size
    FldQty % celsiz = celsiz

    ! Width, height, depth
    allocate(FldQty % whd(size(whd)))
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

  subroutine InitFlowQties3D(FldSlv, xi, values)
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
  impure elemental subroutine InitFlow2D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec2D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1
    real(dp)       , intent(in)    :: value
    integer                        :: ixi(4), AuxWHD(2), x, y, AuxY

    ! (x0,y1) ----------- (x1,y1)
    !         |         |
    !         |         |
    !         |         |
    ! (x0,y0) ----------- (x1,y0)

    ! Set up the simulation bounds.
    bounds: do x = 1, 2
      y = 2*x
      ixi(y-1: y) = aint( xi % x(y-1: y)  / FldQty % celsiz - FldQty % ost(x) )
    end do bounds

    AuxWHD = [max(ixi(1), 0), min(ixi(2), FldQty % whd(1)) - 1]

    ly: do y = max(ixi(3), 0), min(ixi(4), FldQty % whd(2)) - 1
      ! Loop through x.
      ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
      AuxY = FldQty % whd(1) * y
      FldQty % old(AuxWHD(1) + AuxY : AuxWHD(2) + AuxY ) = value
    end do ly

  end subroutine InitFlow2D
!==============================================================================!

!==============================================================================!
  impure elemental subroutine InitFlow3D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec3D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1, xi%x(5) = z0, xi%x(6) = z1
    real(dp)       , intent(in)    :: value
    integer                        :: ixi(6), AuxWHD(4), x, y, z, AuxZ, AuxYZ!, idx

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

    AuxWHD = [max(ixi(3), 0), min(ixi(4), FldQty % whd(2)) - 1,max(ixi(1), 0), min(ixi(2), FldQty % whd(1)) - 1]

    lz: do z = max(ixi(5), 0), min(ixi(6), FldQty % whd(3)) - 1
      ! Loop through y.
      ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
      AuxZ = FldQty % whd(1) * FldQty % whd(2)*z
      ly: do y = AuxWHD(1), AuxWHD(2)
        ! Loop through x.
        ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
        AuxYZ = FldQty % whd(1) * y + AuxZ
        FldQty % old(AuxWHD(3) + AuxYZ : AuxWHD(4) + AuxYZ) = value
      end do ly
    end do lz
  end subroutine InitFlow3D
!==============================================================================!

!==============================================================================!
  subroutine BuildRhs2D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale     ! Scaling factor.
    integer                        :: Auxidx(3), AuxWHD, idx, y, x ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz
    ! Auxiliary X limit.
    AuxWHD = FldSlv % whd(1) - 1
    idx = 0
    ly: do y = 0, FldSlv % whd(2) - 1
      Auxidx(1) = FldSlv % u(1) % whd(1) * y
      Auxidx(2) = FldSlv % u(2) % whd(1) * (y+1)
      Auxidx(3) = FldSlv % u(2) % whd(1) * y
      lx: do x = 0, AuxWHD
        FldSlv % prhs(idx) = -scale * (                                      &
                                        FldSlv % u(1) % old(x+1 + Auxidx(1)) &
                                      - FldSlv % u(1) % old(x   + Auxidx(1)) &
                                      + FldSlv % u(2) % old(x   + Auxidx(2)) &
                                      - FldSlv % u(2) % old(x   + Auxidx(3)) &
                                      )
      idx = idx + 1
      end do lx
    end do ly
  end subroutine BuildRhs2D
  !==============================================================================!

  !==============================================================================!
  subroutine BuildRhs3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale
    integer                        :: Auxidx(8), AuxWHD(2), idx, z, y, x ! Counters and index.

    ! Calculate the scale.
    scale =  1._dp/FldSlv % celsiz
    ! Auxiliary limits.
    AuxWHD = FldSlv % whd - 1
    idx = 0
    lz: do z = 0, FldSlv % whd(3) - 1
      Auxidx(1) = FldSlv % u(1) % whd(1) * FldSlv % u(1) % whd(2) * z
      Auxidx(2) = FldSlv % u(2) % whd(1) * FldSlv % u(2) % whd(2) * z
      Auxidx(3) = FldSlv % u(3) % whd(1) * FldSlv % u(3) % whd(2) * (z+1)
      Auxidx(4) = FldSlv % u(3) % whd(1) * FldSlv % u(3) % whd(2) * z
      ly: do y = 0, AuxWHD(2)
        Auxidx(5) = FldSlv % u(1) % whd(1) * y
        Auxidx(6) = FldSlv % u(2) % whd(1) * (y+1)
        Auxidx(7) = FldSlv % u(3) % whd(1) * y
        Auxidx(8) = FldSlv % u(2) % whd(1) * y
        lx: do x = 0, AuxWHD(1)
          FldSlv % prhs(idx) = -scale * (                                                   &
                                          FldSlv % u(1) % old(x+1 + Auxidx(1) + Auxidx(5) ) &
                                        - FldSlv % u(1) % old(x   + Auxidx(1) + Auxidx(5) ) &
                                        + FldSlv % u(2) % old(x   + Auxidx(6) + Auxidx(2) ) &
                                        - FldSlv % u(2) % old(x   + Auxidx(2) + Auxidx(8) ) &
                                        + FldSlv % u(3) % old(x   + Auxidx(7) + Auxidx(3) ) &
                                        - FldSlv % u(3) % old(x   + Auxidx(4) + Auxidx(7) ) &
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
    real(dp)       , intent(in)    :: DErr   ! Desired error threshold.
    real(dp)                       :: scale, Err, err2 ! Scaling factor, iteration error.
    real(dp)                       :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                       :: NewP ! New pressure value.
    integer                        :: AuxWHD(2)
    integer        , intent(in)    :: MaxIt  ! Maximum number of iterations.
    integer                        :: iter, idx, x, y ! Iteration, counters, index

    ! Set scale.
    scale  = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
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
  subroutine GSPSlv3D(FldSlv, MaxIt, DErr)
    ! Gauss-Seidel pressure solve.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)     , intent(in)      :: DErr   ! Desired error threshold.
    real(dp)                       :: scale, Err ! Scaling factor, iteration error.
    real(dp)                       :: Diag, OffDiag ! Diagonal and offdiagonal elements of the matrix.
    real(dp)                       :: NewP ! New pressure value.
    integer                        :: AuxWHD(3)
    integer      , intent(in)      :: MaxIt  ! Maximum number of iterations.
    integer                        :: iter, idx, z, y, x ! Iteration, counters, index


    ! Set scale.
    scale  = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz * FldSlv % celsiz)
    ! Set range for loops.
    AuxWHD = FldSlv % whd - 1

    ! ITERation Loop.
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
  function Lerp1D(FldQty, yi, x)
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
  impure elemental function Lerp2D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: VtxVal(4) ! VALues at cell VERtices. VtxVal(1) = x00, VtxVal(2) = x10, VtxVal(3) = x01, VtxVal(4) = x11
    type(Vec2D)    , intent(in) :: xi ! xi%x(1) = x, xi%x(2) = y
    real(dp)                    :: x(2), OneMinusX(2)
    real(dp)                    :: Lerp2D
    integer, parameter          :: ivtc(8) = [0,1,0,1,0,0,1,1] ! Index array of the cell vertex displacements.
    integer                     :: ix(2), i

    ! x01   -----------   x11
    !       |         |
    !       |         |
    !       |         |
    ! x00   -----------   x10

    ! Clamp coordinates to the boundaries.
    x = min( max(xi % x - FldQty % ost(1:2), 0._dp), real(FldQty % whd - 1.001_dp, dp) )

    ! Prepare interpolation values.
    ix = aint( x )

    ! We make sure 0 < x < 1
    x  = x - ix

    ! Values at vertices.
    do i = 1, 4
      VtxVal(i) = FldQty % old(ix(1) + ivtc(i) + FldQty % whd(1)*(ix(2) + ivtc(i+4)))
    end do

    OneMinusX = 1._dp-x
    Lerp2D = VtxVal(1)*(OneMinusX(1))*(OneMinusX(2)) + &
             VtxVal(2)*x(1)*(OneMinusX(2))           + &
             VtxVal(3)*x(2)*(OneMinusX(1))           + &
             VtxVal(4)*product(x)
  end function Lerp2D
!==============================================================================!

!==============================================================================!
  impure elemental function Lerp3D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: VtxVal(8) ! Values at cell vertices. VtxVal(1) = x000, VtxVal(2) = x100, VtxVal(3) = x010, VtxVal(4) = x110, VtxVal(5) = x001, VtxVal(6) = x101, VtxVal(7) = x011, VtxVal(8) = x111
    type(Vec3D),     intent(in) :: xi! xi%x(1) = x, xi%x(2) = y, xi%x(3) = z
    real(dp)                    :: x(3), OneMinusX(3)
    real(dp)                    :: Lerp3D
    integer, parameter          :: ivtc(24) = [0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1] ! Index array of the cell vertex displacements.
    integer                     :: ix(3), i

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
    x = min( max(xi % x - FldQty % ost(1:3), 0._dp), real(FldQty % whd - 1.001_dp, dp) )

    ! Prepare interpolation values.
    ix = aint( x )
    ! We make sure 0 < x < 1
    x  = x - ix

    ! Values at vertices.
    do i = 1, 8
      VtxVal(1) = FldQty % old(ix(1)+ivtc(i)+FldQty%whd(1)*(ix(2)+ivtc(i+8))+FldQty%whd(1)*FldQty%whd(2)*(ix(3)+ivtc(i+16)))
    end do

    OneMinusX = 1._dp-x
    Lerp3D = VtxVal(1)*product(OneMinusX)             + &
             VtxVal(2)*x(1)*product(OneMinusX(2:3))   + &
             VtxVal(3)*x(2)*product(OneMinusX(1:3:2)) + &
             VtxVal(4)*product(x(1:2))*OneMinusX(3)   + &
             VtxVal(5)*x(3)*product(OneMinusX(1:2))   + &
             VtxVal(6)*product(x(1:3:2))*OneMinusX(2) + &
             VtxVal(7)*product(x(2:3))*OneMinusX(1)   + &
             VtxVal(8)*product(x)

    !Lerp3D = FldQty % Lerp(                            &
    !            [                                      &
    !              ! Front face x000, x100, x010, x110
    !              FldQty % Lerp(                       &
    !                [                                  &
    !                  ! Lower edge: x000, x100.
    !                  FldQty % Lerp(VtxVal(1:2),x(1)), &
    !                  ! Upper edge: x010, x110.
    !                  FldQty % Lerp(VtxVal(3:4),x(1))  &
    !                ]                                , &
    !                x(2)                               &
    !                           ),                      &
    !              ! Back face x001, x101, x011, x111
    !              FldQty % Lerp(                       &
    !                [                                  &
    !                  ! Lower edge: x001, x101.
    !                  FldQty % Lerp(VtxVal(5:6),x(1)), &
    !                  ! Upper edge: x011, x111.
    !                  FldQty % Lerp(VtxVal(7:8),x(1))  &
    !                ]                                , &
    !                x(2)                               &
    !                          )                        &
    !            ],                                     &
    !            x(3)                                   &
    !                      )
  end function Lerp3D
!==============================================================================!

!==============================================================================!
  subroutine MaxTimestep2D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, SRVel ! Eucledian norm of the velocity, Max(Vel)
    real(dp)                       :: cy, cx ! Center i, center j
    integer                        :: AuxWHD, y, x   ! Counters

    MaxVel = 0._dp
    AuxWHD = FldSlv % whd(1) - 1
    ly: do y = 0, FldSlv % whd(2) - 1
      ! Half way to the next y-index.
      cy = y + 0.5_dp
      lx: do x = 0, AuxWHD
        ! Half way to the next x-index.
        cx = x + 0.5_dp
        ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
        SRVel = sum( FldSlv % u % Lerp(Vec2D([cx, cy])) ** 2 )
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
!==============================================================================!

  subroutine MaxTimestep3D(FldSlv, t)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)       , intent(in)    :: t ! Time
    real(dp)                       :: MaxVel, SRVel ! Eucledian norm of the velocity, Max(Vel).
    real(dp)                       :: cz, cy, cx ! Center i, center j, center k
    integer                        :: AuxWHD(2), z, y, x    ! Counters

    MaxVel = 0._dp
    AuxWHD = FldSlv % whd - 1
    lz: do z = 0, FldSlv % whd(3) - 1
      cz = z + 0.5_dp
      ly: do y = 0, AuxWHD(2)
        ! Half way to the next y-index.
        cy = y + 0.5_dp
        lx: do x = 0, AuxWHD(1)
          ! Half way to the next x-index.
          cx = x + 0.5_dp
          ! Calculate the Root-mean-squared * sqrt(n). Dividing by sqrt(n) is unecessary because n >= 1.
          SRVel = sum( FldSlv % u % Lerp(Vec3D([cx, cy, cz])) ** 2 )
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

  subroutine FEulerVel2D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(FluidQty) , intent(in)  :: FldVel(2)
    class(FluidQty), intent(in)  :: FldQty
    type(Vec2D)    , intent(in)  :: xi
    real(dp)       , intent(out) :: xf(2)
    real(dp)                     :: vel(2) ! Velocities.
    real(dp)       , intent(in)  :: h

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    xf = xi % x - vel * h
  end subroutine FEulerVel2D

  subroutine FEulerVel3D(FldQty, FldVel, xi, xf, h)
    implicit none
    type(FluidQty) , intent(in)  :: FldVel(3)
    class(FluidQty), intent(in)  :: FldQty
    type(Vec3D)    , intent(in)  :: xi
    real(dp)       , intent(out) :: xf(3)
    real(dp)                     :: vel(3) ! Velocities.
    real(dp)       , intent(in)  :: h

    vel = FldVel % Lerp(xi) / FldQty % celsiz

    xf = xi % x - vel * h
  end subroutine FEulerVel3D

  subroutine Advect2D(FldQty, FldVel, h)
    implicit none
    type(FluidQty2D), intent(in)    :: FldVel
    class(FluidQty) , intent(inout) :: FldQty
    real(dp)                        :: xi(2), xf(2)
    real(dp)        , intent(in)    :: h
    integer                         :: AuxWHD, idx, y, x

    AuxWHD = FldQty % whd(1) - 1
    idx = 0
    ! Loop through y.
    ly: do y = 0, FldQty % whd(2) - 1
      xi(2) = y + FldQty % ost(2)
      ! Loop through x.
      lx: do x = 0, AuxWHD
        xi(1) = x + FldQty % ost(1)
        ! Integrate FldQty.
        call FldQty % FEuler(FldVel % FldQty, Vec2D(xi), xf, h)
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
    integer                         :: AuxWHD(2), idx, z, y, x

    AuxWHD = FldQty % whd - 1
    idx = 0
    ! Loop through z.
    lz: do z = 0, FldQty % whd(3) - 1
      xi(3) = z + FldQty % ost(3)
      ! Loop through y.
      ly: do y = 0, AuxWHD(2)
        xi(2) = y + FldQty % ost(2)
        ! Loop through x.
        lx: do x = 0, AuxWHD(1)
          xi(1) = x + FldQty % ost(1)
          ! Integrate FldQty.
          call FldQty % FEuler(FldVel % FldQty, Vec3D(xi), xf, h)
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
    integer                        :: Auxidx(7), idx, y, x

    scale = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz)

    AuxWHD = FldSlv % whd - 1

    ! Applies pressure and updates the velocity field.
    idx = 0
    ly: do y = 0, AuxWHD(2)
      Auxidx(1) = FldSlv % u(1) % whd(1) * y
      Auxidx(2) = FldSlv % u(2) % whd(1) * y
      Auxidx(3) = FldSlv % u(2) % whd(1) * (y+1) ! 5 <-> 3
      lx: do x = 0, AuxWHD(1)
        Auxidx(4:5) = Auxidx(1:2) + x
        Auxidx(6)   = Auxidx(1) + x+1
        Auxidx(7)   = Auxidx(3) + x
        sxp = scale * FldSlv % p(idx)
        FldSlv % u(1) % old(Auxidx(4)) = FldSlv % u(1) % old(Auxidx(4)) - sxp
        FldSlv % u(2) % old(Auxidx(5)) = FldSlv % u(2) % old(Auxidx(5)) - sxp
        FldSlv % u(1) % old(Auxidx(6)) = FldSlv % u(1) % old(Auxidx(6)) + sxp
        FldSlv % u(2) % old(Auxidx(7)) = FldSlv % u(2) % old(Auxidx(7)) + sxp
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
    vb: do y = 0, AuxWHD(2)
      Auxidx(1) = FldSlv % u(1) % whd(1) * y
      ! Lower x boundary.
      FldSlv % u(1) % old(Auxidx(1)) = 0._dp
      ! Upper x boundary
      FldSlv % u(1) % old(FldSlv % whd(1) + Auxidx(1)) = 0._dp
    end do vb

    ! Horizontal boundaries.
    Auxidx(1) = FldSlv % u(2) % whd(1) * FldSlv % whd(2)
    hb: do x = 0, AuxWHD(1)
      ! Lower y boundary.
      FldSlv % u(2) % old(x) = 0._dp
      ! Upper y boundary
      FldSlv % u(2) % old(x + Auxidx(1)) = 0._dp
    end do hb

  end subroutine ApplyPressure2D

  subroutine ApplyPressure3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale, sxp ! Scale, scale x pressure.
    integer                        :: AuxWHD(3)
    integer                        :: Auxidx(15), idx, z, y, x, yp1, xp1

    scale = FldSlv % TimeStep(1) / (FldSlv % srho * FldSlv % celsiz)

    AuxWHD = FldSlv % whd - 1

    ! Applies pressure and updates the velocity field.
    idx = 0
    lz: do z = 0, AuxWHD(3)
      Auxidx(1) = FldSlv % u(1) % whd(1) * FldSlv % u(1) % whd(2) * z
      Auxidx(2) = FldSlv % u(2) % whd(1) * FldSlv % u(2) % whd(2) * z
      Auxidx(3) = FldSlv % u(3) % whd(1) * FldSlv % u(3) % whd(2) * z
      Auxidx(4) = FldSlv % u(3) % whd(1) * FldSlv % u(3) % whd(2) * (z+1)
      ly: do y = 0, AuxWHD(2)
        Auxidx(5) = FldSlv % u(1) % whd(1) * y + Auxidx(1)
        Auxidx(6) = FldSlv % u(2) % whd(1) * y + Auxidx(2)
        Auxidx(7) = FldSlv % u(3) % whd(1) * y + Auxidx(3)
        Auxidx(8) = FldSlv % u(2) % whd(1) * (y+1) + Auxidx(2)
        Auxidx(9) = FldSlv % u(3) % whd(1) * y + Auxidx(4)
        lx: do x = 0, AuxWHD(1)
          Auxidx(10) = x + Auxidx(5)
          Auxidx(11) = x + Auxidx(6)
          Auxidx(12) = x + Auxidx(7)
          Auxidx(13) = x+1 + Auxidx(5)
          Auxidx(14) = x + Auxidx(8)
          Auxidx(15) = x + Auxidx(9)
          xp1 = x + 1
          sxp = scale * FldSlv % p(idx)
          ! Update velocity field.
          FldSlv % u(1) % old(Auxidx(10)) = FldSlv % u(1) % old(Auxidx(10)) - sxp
          FldSlv % u(2) % old(Auxidx(11)) = FldSlv % u(2) % old(Auxidx(11)) - sxp
          FldSlv % u(3) % old(Auxidx(12)) = FldSlv % u(3) % old(Auxidx(12)) - sxp
          FldSlv % u(1) % old(Auxidx(13)) = FldSlv % u(1) % old(Auxidx(13)) + sxp
          FldSlv % u(2) % old(Auxidx(14)) = FldSlv % u(2) % old(Auxidx(14)) + sxp
          FldSlv % u(3) % old(Auxidx(15)) = FldSlv % u(3) % old(Auxidx(15)) + sxp
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
    lz2: do z = 0, AuxWHD(3)
      Auxidx(1) = FldSlv % u(1) % whd(1) * FldSlv % u(1) % whd(2) * z
      ! Loop through y.
      ly3: do y = 0, AuxWHD(2)
        Auxidx(2) = FldSlv % u(1) % whd(1) * y + Auxidx(1)
        ! Left YZ plane (ix = 0)
        FldSlv % u(1) % old(Auxidx(1)) = 0._dp
        ! Right YZ plane (ix = whd(1))
        FldSlv % u(1) % old(FldSlv % whd(1) + Auxidx(1)) = 0._dp
      end do ly3
    end do lz2

    ! XY Plane Boundaries
    ! Loop through y.
    Auxidx(2) = FldSlv % u(3) % whd(1) * FldSlv % u(3) % whd(2) * FldSlv % whd(3)
    ly2: do y = 0, AuxWHD(2)
      Auxidx(1) = FldSlv % u(3) % whd(1) * y
      ! Loop through x.
      lx2: do x = 0, AuxWHD(1)
        ! Front XY plane (iz = 0)
        FldSlv % u(3) % old(x + Auxidx(1)) = 0._dp
        ! Back XY plane (iz = whd(3))
        FldSlv % u(3) % old(x + Auxidx(1) + Auxidx(2)) = 0._dp
      end do lx2
    end do ly2

    ! XZ Plane Boundaries.
    ! Loop through z.
    Auxidx(2) = FldSlv % u(2) % whd(1) * FldSlv % whd(2)
    lz3: do z = 0, AuxWHD(3)
      Auxidx(1) = FldSlv % u(2) % whd(1) * FldSlv % u(2) % whd(2) * z
      ! Loop through x.
      lx3: do x = 0, AuxWHD(1)
        ! Bottom XZ plane (iy = 1)
        FldSlv % u(2) % old(x + Auxidx(1)) = 0._dp
        ! Top XZ plane (iy = whd(2))
        FldSlv % u(2) % old(x + Auxidx(1) + Auxidx(2)) = 0.
        !call FldSlv % u(2) % WriteOldVal(x, FldSlv % whd(2), z, 0._dp)
      end do lx3
    end do lz3
  end subroutine ApplyPressure3D

  function ReadArray2D(FldSlv, array, i, j)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    real(dp)                    :: ReadArray2D
    integer,         intent(in) :: i, j

    ReadArray2D = array(i + FldSlv % whd(1)*j)
  end function ReadArray2D

  function ReadArray3D(FldSlv, array, i, j, k)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(in) :: array(0:)
    real(dp)                    :: ReadArray3D
    integer,         intent(in) :: i, j, k

    ReadArray3D = array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*FldSlv % whd(2)*k)
  end function ReadArray3D

  subroutine WriteArray2D(FldSlv, array, i, j, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    real(dp),        intent(in) :: value
    integer,         intent(in) :: i, j

    array(i + FldSlv % whd(1)*j) = value
  end subroutine WriteArray2D

  subroutine WriteArray3D(FldSlv, array, i, j, k, value)
    implicit none
    class(FluidSlv), intent(in) :: FldSlv
    real(dp),        intent(inout) :: array(0:)
    real(dp),        intent(in) :: value
    integer,         intent(in) :: i, j, k

    array(i + FldSlv % whd(1)*j + FldSlv % whd(1)*FldSlv % whd(2)*k) = value
  end subroutine WriteArray3D

  impure elemental function ReadVal2D(FldQty, i, j)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal2D
    integer,         intent(in) :: i, j

    ReadVal2D = FldQty % old(i + FldQty % whd(1)*j)
  end function ReadVal2D

  impure elemental function ReadVal3D(FldQty, i, j, k)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal3D
    integer,         intent(in) :: i, j, k

    ReadVal3D = FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1)*FldQty % whd(2)*k)
  end function ReadVal3D

  impure elemental subroutine WriteNewVal2D(FldQty, i, j, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j

    FldQty % new(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteNewVal2D

  impure elemental subroutine WriteNewVal3D(FldQty, i, j, k, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j, k

    FldQty % new(i + FldQty % whd(1)*j + FldQty % whd(1) * FldQty % whd(2) * k) = NewVal
  end subroutine WriteNewVal3D

  impure elemental subroutine WriteOldVal2D(FldQty, i, j, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j

    FldQty % old(i + FldQty % whd(1)*j) = NewVal
  end subroutine WriteOldVal2D

  impure elemental subroutine WriteOldVal3D(FldQty, i, j, k, NewVal)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    real(dp),        intent(in)    :: NewVal
    integer,         intent(in)    :: i, j, k

    FldQty % old(i + FldQty % whd(1)*j + FldQty % whd(1)*FldQty % whd(2)*k) = NewVal
  end subroutine WriteOldVal3D

  impure elemental subroutine UpdateVals(FldQty)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    FldQty % old = FldQty % new
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
