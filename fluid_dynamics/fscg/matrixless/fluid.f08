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
    procedure :: AddFlow2D, AddFlow3D
    generic   :: AddFlow => AddFlow2D, AddFlow3D
    ! Read values.
    procedure :: ReadVal2D, ReadVal3D
    generic   :: ReadVal => ReadVal2D, ReadVal3D
    ! Write values.
    procedure :: WriteNewVal2D, WriteNewVal3D
    generic   :: WriteNewVal => WriteNewVal2D, WriteNewVal3D
    procedure :: WriteOldVal2D, WriteOldVal3D
    generic   :: WriteOldVal => WriteOldVal2D, WriteOldVal3D
    ! Update values.
    procedure :: UpdateVals
    ! Linear interpolate.
    procedure :: Lerp1D, Lerp2D, Lerp3D
    generic   :: Lerp => Lerp1D, Lerp2D, Lerp3D
    ! Integrate.
    procedure :: FEulerVel2D, FEulerVel3D
    generic   :: FEuler => FEulerVel2D, FEulerVel3D
    ! Advect.
    procedure :: Advect2D, Advect3D
    generic   :: Advect => Advect2D, Advect3D
  end type FluidQty

  type FluidSlv
    character(len=:), dimension(:), allocatable :: names ! Names of variables.
    type(FluidQty) , allocatable :: u(:)   ! Velocities
    type(FluidQty)               :: rho    ! Density array
    integer        , allocatable :: whd(:) ! Width, height, depth
    real(dp)                     :: celsiz ! Cell size.
    real(dp)                     :: irho   ! Initial density.
    real(dp)       , allocatable :: prhs(:)! Right hand side of pressure equation.
    real(dp)       , allocatable :: p(:)   ! Pressure solution.
  contains
    ! Initialise type.
    procedure :: Init => InitFluid
    ! Build right hand side.
    procedure :: BuildRhs, BuildRhs2D, BuildRhs3D
  end type FluidSlv

  type Vec2D
    real(dp) :: x(2)
  end type Vec2D

  type FluidQty2D
    type(FluidQty) :: FldQty(2)
  end type FluidQty2D

  type Vec2D2
    real(dp) :: x(4)
  end type Vec2D2

  type Vec3D
    real(dp) :: x(3)
  end type Vec3D

  type FluidQty3D
    type(FluidQty) :: FldQty(3)
  end type FluidQty3D

  type Vec3D2
    real(dp) :: x(6)
  end type Vec3D2

contains


  pure function Lerp1D(FldQty, yi, x)
    ! Linear intERPolate in 1D.
    ! 1D interpolation between yi(1) < y < yi(4) for 0 < x < 1
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)       , intent(in) :: yi(2), x
    real(dp)                    :: lerp1d

    Lerp1D = yi(1)*(1.0 - x) + yi(2)*x

  end function Lerp1D

  elemental function Lerp2D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    type(Vec2D)    , intent(in) :: xi ! xi%x(1) = x, xi%x(2) = y
    real(dp)                    :: Lerp2D
    real(dp)                    :: x(2)
    integer                     :: ix(2)
    real(dp)                    :: VtxVal(4) ! VALues at cell VERtices. VtxVal(1) = x00, VtxVal(2) = x10, VtxVal(3) = x01, VtxVal(4) = x11
    integer, parameter          :: ivtc(8) = [0,1,0,1,0,0,1,1] + 1 ! Index array of the cell vertex displacements.
    ! x01   -----------   x11
    !       |         |
    !       |         |
    !       |         |
    ! x00   -----------   x10

    ! Clamp coordinates to the boundaries.
    x = min( max(xi % x - FldQty % ost(1:2), 0.), float(FldQty % whd(1:2)) )

    ! Prepare interpolation values.
    ix = floor( x )
    ! We make sure 0 < x < 1
    x  = x - ix

    ! Values at vertices.
    VtxVal(1:4) = FldQty % ReadVal(ix(1) + ivtc(1:4), ix(2) + ivtc(5:8))

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

  elemental function Lerp3D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    type(Vec3D),     intent(in) :: xi! xi%x(1) = x, xi%x(2) = y, xi%x(3) = z
    real(dp)                    :: Lerp3D
    real(dp)                    :: x(3)
    integer                     :: ix(3)
    real(dp)                    :: VtxVal(8) ! Values at cell vertices. VtxVal(1) = x000, VtxVal(2) = x100, VtxVal(3) = x010, VtxVal(4) = x110, VtxVal(5) = x001, VtxVal(6) = x101, VtxVal(7) = x011, VtxVal(8) = x111
    integer, parameter          :: ivtc(24) = [0,1,0,1,0,1,0,1,0,0,1,1,0,0,1,1,0,0,0,0,1,1,1,1] + 1 ! Index array of the cell vertex displacements.

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
    x = min( max(xi % x - FldQty % ost(1:3), 0.), float(FldQty % whd(1:3)) )

    ! Prepare interpolation values.
    ix = floor( x )
    ! We make sure 0 < x < 1
    x  = x - ix

    VtxVal(1:8) = FldQty % ReadVal(ix(1) + ivtc(1:8), ix(2) + ivtc(9:16), ix(3) + ivtc(17:24))

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
    integer                         :: i, j

    ! Loop through y.
    ly: do j = 1, FldQty % whd(2)
      ! Set value of y.
      x(2) = j + FldQty % ost(2)
      ! Loop through x.
      lx: do i = 1, FldQty % whd(1)
        x(1) = i + FldQty % ost(1)
        ! Integrate FldQty.
        call FldQty % FEuler([FldVel % FldQty(1), FldVel % FldQty(2)], Vec2D(x), x, h)
        ! Interpolate component from grid and write new value.
        call FldQty % WriteNewVal(i, j, FldQty % Lerp(Vec2D(x)))
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
    lz: do k = 1, FldQty % whd(3)
      ! Set value of z.
      x(3) = k + FldQty % ost(3)
      ! Loop through y.
      ly: do j = 1, FldQty % whd(2)
        ! Set value of y.
        x(2) = j + FldQty % ost(2)
        ! Loop through x.
        lx: do i = 1, FldQty % whd(1)
          x(1) = i + FldQty % ost(1)
          ! Integrate FldQty.
          call FldQty % FEuler([FldVel % FldQty(1), FldVel % FldQty(2), FldVel % FldQty(3)], Vec3D(x), x, h)
          ! Interpolate component from grid and save the new value.
          call FldQty % WriteNewVal(i, j, k, FldQty % Lerp(Vec3D(x)))
        end do lx
      end do ly
    end do lz
  end subroutine Advect3D

  subroutine AddFlow2D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec2D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1
    real(dp)       , intent(in)    :: value
    integer                        :: i, j, ixi(4)

    ! (x0,y1) ----------- (x1,y1)
    !         |         |
    !         |         |
    !         |         |
    ! (x0,y0) ----------- (x1,y0)

    ! Set up the simulation boundaries.
    bounds: do i = 1, 2
      j = 2*i
      ixi(j-1 : j) = floor( xi%x(j-1 : j) / FldQty % celsiz - FldQty % ost(i) ) + 1
    end do bounds

    ! Loop through y.
    ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
    ly: do j = max(ixi(3), 1), min( ixi(4), FldQty % whd(2) )
      ! Loop through x.
      ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
      lx: do i = max(ixi(1), 1), min( ixi(2), FldQty % whd(2) )
        ! Write value to old value.
        if ( abs(FldQty % ReadVal(i, j)) < abs(value) ) call FldQty % WriteOldVal(i,j,value)
      end do lx
    end do ly

  end subroutine AddFlow2D

  subroutine AddFlow3D(FldQty, xi, value)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    type(Vec3D2)   , intent(in)    :: xi ! xi%x(1) = x0, xi%x(2) = x1, xi%x(3) = y0, xi%x(4) = y1, xi%x(5) = z0, xi%x(6) = z1
    real(dp)       , intent(in)    :: value
    integer                        :: i, j, k, ixi(6)

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
    bounds: do i = 1, 3
      j = 2*i
      ixi(j-1 : j) = floor( xi%x(j-1 : j) / FldQty % celsiz - FldQty % ost(i) )
    end do bounds

    ! Loop through z.
    ! Start from the maximum value between the desired z0 and 1 (for the cell index) up to the minimum value between the desired z1 and the simulation depth.
    lz: do k = max(ixi(5), 1), min( ixi(6), FldQty % whd(3) )
      ! Loop through y.
      ! Start from the maximum value between the desired y0 and 1 (for the cell index) up to the minimum value between the desired y1 and the simulation height.
      ly: do j = max(ixi(3), 1), min( ixi(4), FldQty % whd(2) )
        ! Loop through x.
        ! Start from the maximum value between the desired x0 and 1 (for the cell index) up to the minimum value between the desired x1 and the simulation width.
        lx: do i = max(ixi(1), 1), min( ixi(2), FldQty % whd(2) )
          ! Write value to old value.
          if ( abs(FldQty % ReadVal(i, j, k)) < abs(value) ) call FldQty % WriteOldVal(i, j, k, value)
        end do lx
      end do ly
    end do lz
  end subroutine AddFlow3D

  subroutine BuildRhs(FldSlv)
    ! Build Rhs wrapper subroutine.
    ! The actual subroutines build the pressure solution via the negative of the divergence as a forward finite difference (truly horrible). This is not the best way to go about doing this, but it will be fixed in the future.
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    integer, allocatable, save     :: ndim ! Number of dimensions.

    if (.not. allocated(ndim)) then
      allocate(ndim, source = size(FldSlv % u))
    end if

    if (ndim == 2) then
      call FldSlv % BuildRhs2D
    else
      call FldSlv % BuildRhs3D
    end if
  end subroutine BuildRhs

  subroutine BuildRhs2D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale
    integer                        :: i, j, idx

    ! Calculate the scale.
    scale =  1./FldSlv % celsiz

    idx = 1
    ly: do j = 1, FldSlv % whd(2)
      lx: do i = 1, FldSlv % whd(1)
        FldSlv % prhs(idx) = -scale * ( &
                                         FldSlv % u(1) % ReadVal(i+1, j  ) &
                                       + FldSlv % u(2) % ReadVal(i  , j+1) &
                                       -sum(FldSlv % u % ReadVal(i  , j  ))&
                                      )
        idx = idx + 1
      end do lx
    end do ly
  end subroutine BuildRhs2D

  subroutine BuildRhs3D(FldSlv)
    implicit none
    class(FluidSlv), intent(inout) :: FldSlv
    real(dp)                       :: scale
    integer                        :: i, j, k, idx

    ! Calculate the scale.
    scale =  1./FldSlv % celsiz

    idx = 1
    lz: do k = 1, FldSlv % whd(3)
      ly: do j = 1, FldSlv % whd(2)
        lx: do i = 1, FldSlv % whd(1)
          FldSlv % prhs(idx) = -scale * ( &
                                          FldSlv % u(1) % ReadVal(i+1, j  , k  )  &
                                        + FldSlv % u(2) % ReadVal(i  , j+1, k  )  &
                                        + FldSlv % u(3) % ReadVal(i  , j  , k+1)  &
                                        -sum(FldSlv % u % ReadVal(i  , j  , k  )) &
                                        )
          idx = idx + 1
        end do lx
      end do ly
    end do lz
  end subroutine BuildRhs3D

  subroutine InitFluid(FldSlv, names, whd, scale)
    implicit none
    character(len=*),  dimension(:), intent(in) :: names
    integer           , intent(inout) :: whd(:) ! Width (i), height (j), depth (k)
    real(dp), optional, intent(in) :: scale
    class(FluidSlv),   intent(out) :: FldSlv ! Fluid
    real(dp), dimension(size(whd)) :: ost    ! Offset
    real(dp), dimension(size(whd)) :: TMPost
    real(dp)                       :: celsiz ! Cell size
    integer                        :: ndim   ! Number of dimensions
    integer                        :: pwhd   ! Product of the dimensions
    integer, dimension(size(whd))  :: TMPwhd ! whd for each different quantity
    integer                        :: i      ! Counter

    ! Number of dimensions.
    ndim = size(whd)
    pwhd = product(whd)

    ! Cell size.
    if(present(scale) .eqv. .true.) then
      celsiz = scale/minval(whd)
    else
      celsiz = 1./minval(whd)
    endif

    ! Initialise density.
    ! Centered in the grid.
    ost    = 0.5
    call FldSlv % rho % init(trim(names(1)(1:len(names(1)))), whd, ost, celsiz)

    ! Initialise velocities.
    ! Staggered on the grid.
    allocate(FldSlv % u(ndim))
    iv: do i = 1, ndim
      TMPwhd    = whd
      TMPwhd(i) = whd(i) + 1
      TMPost    = ost
      TMPost(i) = ost(i) - 0.5
      call FldSlv % u(i) % init(trim(names(i+1)(1:len(names(i+1)))), TMPwhd, TMPost, celsiz)
    end do iv

    ! Allocate whd (simulation bounds).
    allocate (FldSlv % whd(ndim))
    FldSlv % whd = whd

    ! Allocate the right hand side of pressure solve.
    allocate(FldSlv % prhs(pwhd))
    FldSlv % prhs = 0.

    ! Allocate the pressure solution.
    allocate(FldSlv % p(pwhd))
    FldSlv % p = 0.

    FldSlv % celsiz = celsiz

  end subroutine InitFluid

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
    ncells = product(whd)
    allocate(FldQty % old(ncells), FldQty % new(ncells))
    FldQty % old = 0.
    FldQty % new = 0.

    ! Cell size
    FldQty % celsiz = celsiz
  end subroutine InitFluidQty

  elemental function ReadVal2D(FldQty, i, j)
    implicit none
    integer,         intent(in) :: i, j
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal2D
    integer                     :: idx

    idx = i + FldQty % whd(1)*(j-1)
    ReadVal2D = FldQty % old(idx)
  end function ReadVal2D

  elemental function ReadVal3D(FldQty, i, j, k)
    implicit none
    integer,         intent(in) :: i, j, k
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal3D
    integer                     :: idx

    idx = i + FldQty % whd(1)*(j-1) + FldQty % whd(1) * FldQty % whd(2)*(k-1)
    ReadVal3D = FldQty % old(idx)
  end function ReadVal3D

  subroutine WriteNewVal2D(FldQty, i, j, NewVal)
    implicit none
    integer,         intent(in)    :: i, j
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1)
    FldQty % new(idx) = NewVal
  end subroutine WriteNewVal2D

  subroutine WriteNewVal3D(FldQty, i, j, k, NewVal)
    implicit none
    integer,         intent(in)    :: i, j, k
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1) + FldQty % whd(1) * FldQty % whd(2)*(k-1)
    FldQty % new(idx) = NewVal
  end subroutine WriteNewVal3D

  subroutine WriteOldVal2D(FldQty, i, j, NewVal)
    implicit none
    integer,         intent(in)    :: i, j
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1)
    FldQty % old(idx) = NewVal
  end subroutine WriteOldVal2D

  subroutine WriteOldVal3D(FldQty, i, j, k, NewVal)
    implicit none
    integer,         intent(in)    :: i, j, k
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1) + FldQty % whd(1) * FldQty % whd(2)*(k-1)
    FldQty % old(idx) = NewVal
  end subroutine WriteOldVal3D

  subroutine UpdateVals(FldQty)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    FldQty % old = FldQty % new
    FldQty % new = 0.
  end subroutine UpdateVals

end module fluid
