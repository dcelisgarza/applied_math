module fluid
  use nrtype
  type FluidQty
    character(len=:), allocatable :: name ! Variable name.
    real(dp),         allocatable :: old(:) ! Old value.
    real(dp),         allocatable :: new(:) ! New value.
    integer,          allocatable :: whd(:) ! Width, height, depth of the simulation domain.
    real(dp),         allocatable :: ost(:) ! Offset.
    real(dp)                      :: celsiz ! Cell size (area or volume).
  contains
    ! Initialise type.
    procedure :: Init => InitFluidQty
    ! Read values.
    procedure :: ReadVal2D
    procedure :: ReadVal3D
    generic   :: ReadVal => ReadVal2D, ReadVal3D
    ! Write values.
    procedure :: WriteVal2D
    procedure :: WriteVal3D
    generic   :: WriteVal => WriteVal2D, WriteVal3D
    ! Update values.
    procedure :: UpdateVals
  end type FluidQty

contains
  subroutine InitFluidQty(FldQty, name, whd, ost, celsiz)
    implicit none
    character(len=*), intent(in)  :: name   ! Variable name.
    integer,          intent(in)  :: whd(:) ! Width (i), height (j), depth (k)
    real(dp),         intent(in)  :: ost(:) ! Offset
    real(dp),         intent(in)  :: celsiz ! Cell size (area or volume)
    class(FluidQty),  intent(out) :: FldQty ! Fluid Quantity
    integer                       :: ncells ! # of cells = prod(whd)

    ! Check that size(ost) == size(ost).
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

  function ReadVal2D(FldQty, i, j)
    implicit none
    integer,         intent(in) :: i, j
    class(FluidQty), intent(in) :: FldQty
    real(dp)                    :: ReadVal2D
    integer                     :: idx

    idx = i + FldQty % whd(1)*(j-1)
    ReadVal2D = FldQty % old(idx)
  end function ReadVal2D

  function ReadVal3D(FldQty, i, j, k)
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

  function Lerp1D(yi,xi,x)
    ! Linear intERPolate in 1D.
    ! 1D interpolation between yi(1) < y < y(2) for x(1) < x < x(2)
    implicit none
    real(dp), intent(in) :: yi(2), xi(2), x
    real(dp) :: lerp1d

    lerp1d = yi(1)*(xi(2)-x) + yi(2)*(x-xi(1))
  end function Lerp1D


end module fluid
