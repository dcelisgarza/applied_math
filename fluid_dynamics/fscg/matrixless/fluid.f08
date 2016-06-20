! Generalised Fortran version of original c++ code by TunnaBrain (Benedikt Bitterli) https://github.com/tunabrain/incremental-fluids.

module fluid
  use nrtype

  type FluidQty
    character(len=:), allocatable :: name   ! Variable name.
    real(dp)        , allocatable :: old(:) ! Old value.
    real(dp)        , allocatable :: new(:) ! New value.
    real(dp)        , allocatable :: whd(:) ! Width, height, depth of the simulation domain.
    real(dp)        , allocatable :: ost(:) ! Offset.
    real(dp)                      :: celsiz ! Cell size (area or volume).
  contains
    ! Initialise type.
    procedure :: Init => InitFluidQty
    ! Read values.
    procedure :: ReadVal2D, ReadVal3D
    generic   :: ReadVal => ReadVal2D, ReadVal3D
    ! Write values.
    procedure :: WriteVal2D, WriteVal3D
    generic   :: WriteVal => WriteVal2D, WriteVal3D
    ! Update values.
    procedure :: UpdateVals
    ! Linear interpolate.
    procedure :: Lerp1D, Lerp2D, Lerp3D
    generic   :: Lerp => Lerp1D, Lerp2D, Lerp3D
  end type FluidQty

  type FluidSlv
    character(len=:), dimension(:), allocatable :: names ! Names of variables.
    type(FluidQty) , allocatable :: u(:)   ! Velocities
    type(FluidQty)               :: rho    ! Density array
    real(dp)       , allocatable :: whd(:) ! Width, height, depth
    real(dp)                     :: celsiz ! Cell size.
    real(dp)                     :: irho   ! Initial density.
    real(dp)       , allocatable :: prhs(:)! Right hand side of pressure equation.
    real(dp)       , allocatable :: p(:)   ! Pressure solution.
  contains
    ! Initialise type.
    procedure :: Init => InitFluid
  end type FluidSlv

  type Vec2D
    real(dp) :: x(2)
  end type Vec2D

  type Vec3D
    real(dp) :: x(3)
  end type Vec3D

contains
  function Lerp1D(FldQty, xi, x)
    ! Linear intERPolate in 1D.
    ! 1D interpolation between xi(3) < y < xi(4) for x(1) < x < x(2)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    real(dp)       , intent(in) :: xi(4), x
    real(dp)                    :: lerp1d

    Lerp1D = xi(3)*(xi(2)-x) + xi(4)*(x-xi(1))

    print*, "xi(1) = ", xi(1)
    print*, "xi(2) = ", xi(2)
    print*, "xi(3) = ", xi(3)
    print*, "xi(4) = ", xi(4)
    print*, "x = ", x

    print*, "xi(3)*(xi(2)-x) = ", xi(3)*(xi(2)-x)
    print*, "xi(4)*(x-xi(1)) = ", xi(4)*(x-xi(1))
    print*, "Lerp1D = ", Lerp1D
  end function Lerp1D

  function Lerp2D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    type(Vec2D)    , intent(in) :: xi ! xi%x(1) = x, xi%x(2) = y
    real(dp)                    :: Lerp2D
    real(dp)                    :: x(2)
    integer                     :: ix(2)
    real(dp)                    :: cvtc(4) ! Cell vertices. cvtc(1) = x00, cvtc(2) = x10, cvtc(3) = x01, cvtc(4) = x11
    integer, parameter          :: ivtc(8) = [0,1,0,1,0,0,1,1] ! Index array of the cell vertices.
    ! x01   -----------   x11
    !       |         |
    !       |         |
    !       |         |
    ! x00   -----------   x10

    ! Clamp coordinates to the boundaries.
    x = min( max(xi % x - FldQty % ost(1:2), 0.), FldQty % whd(1:2) )

    ! Prepare interpolation values.
    ix = floor( x )
    ! We make sure 0 < x < 1
    x  = x - ix

    cvtc(1:4) = FldQty % ReadVal(ix(1) + ivtc(1:4) + 1, ix(2) + ivtc(5:8) + 1)

    Lerp2D = FldQty % Lerp(                                                 &
                            [                                               &
                              0._dp, 1._dp                                , &
                              FldQty % Lerp([0._dp, 1._dp,cvtc(1:2)],x(1)), &
                              FldQty % Lerp([0._dp, 1._dp,cvtc(3:4)],x(1))  &
                            ]                                            ,  &
                            x(2)                                            &
                           )


  end function Lerp2D

  function Lerp3D(FldQty, xi)
    implicit none
    class(FluidQty), intent(in) :: FldQty
    type(Vec3D),     intent(in) :: xi(3)! xi(1) = x, xi(2) = y, xi(3) = z
    real(dp)                    :: Lerp3D

    ! Make lerp 3D
  end function Lerp3D

  subroutine FEulerVel2D(FldQty, xi, xf, h)
    implicit none
    type(Vec2D), intent(in)           :: xi
    real(dp), intent(in)              :: h
    class(FluidQty), intent(in)       :: FldQty(:)
    real(dp), intent(out)             :: xf(:)
    real(dp), dimension(size(FldQty)) :: vel ! Velocities.
    integer :: i

    ! Loop through Velcoities.
    loop_vel: do i = 1, size(FldQty)
      vel(i) = FldQty(i) % Lerp(xi) / FldQty(i) % celsiz
    end do loop_vel

    xf = xi % x - vel*h

  end subroutine FEulerVel2D

  subroutine InitFluid(FldSlv, names, whd, scale)
    implicit none
    character(len=*),  dimension(:), intent(in) :: names
    integer,           intent(in)  :: whd(:) ! Width (i), height (j), depth (k)
    real(dp), optional, intent(in) :: scale
    class(FluidSlv),   intent(out) :: FldSlv ! Fluid
    real(dp), dimension(size(whd)) :: ost   ! Offset
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

    ! Allocate the right hand side of pressure solve.
    allocate(FldSlv % prhs(pwhd))
    FldSlv % prhs = 0.

    ! Allocate the pressure solution.
    allocate(FldSlv % p(pwhd))
    FldSlv % p = 0.

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

  subroutine WriteVal2D(FldQty, i, j, NewVal)
    implicit none
    integer,         intent(in)    :: i, j
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1)
    FldQty % new(idx) = NewVal
  end subroutine WriteVal2D

  subroutine WriteVal3D(FldQty, i, j, k, NewVal)
    implicit none
    integer,         intent(in)    :: i, j, k
    real(dp),        intent(in)    :: NewVal
    class(FluidQty), intent(inout) :: FldQty
    integer                        :: idx

    idx = i + FldQty % whd(1)*(j-1) + FldQty % whd(1) * FldQty % whd(2)*(k-1)
    FldQty % new(idx) = NewVal
  end subroutine WriteVal3D

  subroutine UpdateVals(FldQty)
    implicit none
    class(FluidQty), intent(inout) :: FldQty
    FldQty % old = FldQty % new
    FldQty % new = 0.
  end subroutine UpdateVals

end module fluid
