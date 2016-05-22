module sfldy
  use ode

  type SFluidQty
    real(dp), allocatable :: w(:) ! Cell width and height and depth.
    real(dp), allocatable :: o(:) ! Offset from top-left of the grid, or top-left-front in 3D.
    real(dp) :: cs                ! Cell size.
  end type SFluidQty

  interface interpolate
    module procedure :: lint
  end interface interpolate
contains


  subroutine sfldyEuler(u, x, h)
    ! Forward Euler for the integration of the velocities.
    implicit none
    class(SFluidQty), intent(in) :: u(:)
    real(dp), intent(in)         :: h ! Fluid vectors (u, v, w), time step
    real(dp), intent(inout)      :: x(:)
    real(dp)                     :: lint
    real(dp)                     :: du(size(u))
    integer                      :: i

    ! Obtain terivatives.
    do i = 1, size(u)
      du(i) = interp(u(i), x) / u(i) % cs
    end do

    ! Integrate x.
    x = x + du * h

  end subroutine sfldyEuler

  function interp(u, x) result(interpolation)
    implicit none
    class(SFluidQty), intent(in)     :: u
    real(dp), intent(in)             :: x(:)
    real(dp)                         :: interpolation
    real(dp), dimension(size(x))     :: xt
    integer,  dimension(size(x))     :: ix
    real(dp), dimension(2**size(x))  :: grid
    real(dp), parameter              :: offset = 1+1d-30

    ! Clamping vertex coordinates.
    xt = min(max(x - u % o, 0._dp), u % w - offset)
    ! Making vertices integers.
    ix = nint(x)
    ! Setting up origin to the top-left or top-left-front corner.
    xt = xt - ix

    ! Find grid positions.
    grid = FindCoord(u, x)

    ! This is for 2D only.
    interpolation = lint( [ lint(grid(1:2), x(1)), lint(grid(3:4), x(1)) ], x(2) )
  end function interp

  function lint(c, x)
    implicit none
    real(dp), intent(in) :: c(:), x ! Coeffcients and x
    real(dp)             :: lint

    lint = ( c(1) - c(2) ) * x - c(1)
  end function lint

  function FindCoord(u, x) result(grid)
    ! Scales the offset to the cell width.
    implicit none
    real(dp), intent(in)         :: x(:)
    class(SFluidQty), intent(in) :: u
    real(dp)                     :: grid(2**size(x)), tmp
    integer                      :: i, j, k, sx, sg
    integer, parameter :: disp(24) = [0,1,0,1,0,1,0,1,0,0,1,1,&
                                      0,0,1,1,0,0,0,0,1,1,1,1]

    sx  = size(x)
    sg  = size(grid)

    ! Loop over Grid Displacements.
    gd: do i = 1, sg
      ! Sum Only Relevant Quantities.
      sorq1: do j = 1, sx
        ! Save the value of x(i) in tmp.
        tmp = x(j) + disp(i + 8*(j-1))
        sorq2: do k = 1, sx
          ! If j = i we need to skip to the next value of j.
          if (j == k) cycle sorq2
          ! Check whether this is correct, do not be an idiot!
          tmp = tmp + ( x(k) + disp(i + 8*(k-1)) ) * u % w(j)
          ! Check whether this is correct, do not be an idiot!
        end do sorq2
        ! Assign tmp to the coordinate i.
        grid(i + 8*(j-1)) = tmp
      end do sorq1
    end do gd
  end function FindCoord

end module sfldy
