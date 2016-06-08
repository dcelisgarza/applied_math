module interpolation
  use nrtype
  ! http://www.paulinternet.nl/?page=bicubic
  ! http://citeseerx.ist.psu.edu/viewdoc/download;jsessionid=5C37DCC31C0A337C28A63330C436A3B5?doi=10.1.1.89.7835&rep=rep1&type=pdf
  ! http://bmia.bmt.tue.nl/people/BRomeny/Courses/8C080/Interpolation.pdf
  ! http://paulbourke.net/miscellaneous/interpolation/
contains
  subroutine pinex(xi,yi,x,y,dy)
    !=========================================================!
    ! Polynomial INterpolation and EXtrapolation (PINEX).     !
    ! Runge-Kutta 4 Gill's method                             !
    ! Daniel Celis Garza 24 Sept. 2015                        !
    !---------------------------------------------------------!
    ! Inputs:                                                 !
    ! x    = independent variable                             !
    ! xi() = independent variable points                      !
    ! yi() = independent variable points                      !
    !---------------------------------------------------------!
    ! Outputs:                                                !
    ! y  = dependent variable                                 !
    ! dy = dependent variable increment/error                 !
    !---------------------------------------------------------!
    ! Locals:                                                 !
    ! c()        = small difference 1                         !
    ! d()        = small difference 2                         !
    ! difx()     = distance between point i and input x       !
    ! aux_difx() = aux var to find distance of x to point i   !
    ! den()      = denominator                                !
    ! idx        = closest point to x                         !
    ! aux_idx    = aux var to find index i                    !
    ! i          = counter                                    !
    ! n          = size(xi)                                   !
    !=========================================================!
    implicit none
    real(dp), intent(in)          :: x, xi(:), yi(:)
    real(dp), intent(out)         :: y, dy
    real(dp), dimension(size(xi)) :: c, d, difx, aux_difx, den
    integer                       :: i, idx, aux_idx(1), n

    ! Making sure we're using equally sized arrays.
    check_size: if(size(xi) /= size(yi)) then
      print*, ' Size of x and y arrays must be equal.'
      return
    endif check_size

    ! The order of the interpolation.
    n = size(xi)

    ! Calculating (x-xi) for all i  [and (xj-x) for all j].
    difx = x-xi

    ! Looking for the index of the entry with the closest value to x from xi.
    aux_difx = abs(x-xi)
    aux_idx  = minloc(aux_difx)
    idx      = aux_idx(1)

    ! Our initial best guess for y.
    y        = yi(idx)

    ! Drop down an index because that's what's required by the relationships of Cs and Ds
    idx      = idx - 1

    ! Initialising relationships of Cs and Ds.
    c = yi
    d = yi

    ! Looping over the columns of Neville's tree.
    nev_loop: do i = 1, n - 1
      den(1:n-i) = difx(1:n-i) - difx(1+i:n)
      if ( any( den(1:n-i) == 0._dp ) ) print*, ' Error: subr: poly_interpol: nev_loop. Division by zero.'
      den(1:n-i) = ( c(2:n-i+1) - d(1:n-i) ) / den(1:n-i)
      d(1:n-i)   = difx(1+i:n) * den(1:n-i)
      c(1:n-i)   = difx(1:n-i) * den(1:n-i)

      ! Calculate error terms.
      err: if (2*idx < n-i) then
        dy  = c(idx+1)
      else err
        dy  = d(idx)
        idx = idx-1
      end if err

      ! Increase the value of y.
      y = y + dy
    end do nev_loop
  end subroutine pinex

  subroutine rinex(xi,yi,x,y,dy)
    ! Polynomial INterpolation and EXtrapolation (PINEX).
    implicit none
    real(dp), intent(in)          :: x, xi(:), yi(:)
    real(dp), intent(out)         :: y, dy
    integer                       :: i, j, idx, aux_idx(1), n
    real(dp), dimension(size(xi)) :: c, d, difx, aux_difx, w, v
    real(dp), parameter           :: tiny = 1.d-25

    ! Making sure we're using equally sized arrays.
    check_size: if(size(xi) /= size(yi)) then
      print*, " Size of x and y/home/krunos/Documents/interpolation/interpol_mod.f08 arrays must be equal."
      return
    end if check_size

    ! The order of the interpolation.
    n = size(xi)

    ! Calculating (x-xi) for all i  [and (xj-x) for all j].
    difx = x-xi

    ! Looking for the index of the entry with the closest value to x from xi.
    aux_difx = abs(x-xi)
    aux_idx  = minloc(aux_difx)
    idx      = aux_idx(1)

    ! Our initial best guess for y.
    y        = yi(idx)

    ! Return when x == xi(idx).
    ret: if (x == xi(idx)) then
      dy = 0._dp
      return
    endif ret

    ! Drop down an index because that's what's required by the relationships of Cs and Ds
    idx      = idx - 1

    ! Initialising relationships of Cs and Ds.
    c = yi
    d = yi + tiny

    ! Looping over the columns of Neville's tree.
    nev_loop: do i = 1, n - 1
      ! Calculating (x - x_{i}) / (x - x_{i+m+1}) * D_{m,i}.
      w( 1: n - i ) = ( x - xi( 1: n - i) ) * d( 1: n - i ) / difx( 1 + i: n )
      ! Calculating (x - x_{i}) / (x - x_{i+m+1}) * D_{m,i} - C_{m,i+1}.
      v( 1: n - i ) = w( 1: n - i ) - c( 2: n - i + 1 )

      ! Making sure we're using equally sized arrays.
      check_underflow: if( any( v( 1: n - 1 ) == 0._dp )  ) then
        print*, " Error: rinex: nev_loop: Denominator = 0, use different points."
        return
      end if check_underflow

      ! Calculate the common factor ( c(m,i+1) - d(m,i) ) / ( (x - x_{i}) / ( x - x_{i+m+1}) * D_{m,i} - C_{m,i+1} ).
      v( 1: n - i ) = ( c(2: n - i + 1) - d( 1: n - i ) ) / v( 1: n - i )
      ! Calcualte D(m+1, i).
      d( 1: n - i ) = c(2: n - i + 1) * v ( 1: n - i )
      ! Calcualte C(m+1, i).
      c( 1: n - i ) = w(1: n - i ) * v ( 1: n - i )

      ! Calculate the interpolated value of y.
      ! Decide by how much we're going to increase y so as to cover Neville's tree the fastest.
      err: if ( 2 * idx < n - i ) then
        dy  = c(idx+1)
      else err
        dy  = d(idx)
        idx = idx-1
      end if err
      y = y + dy
    end do nev_loop
  end subroutine rinex

  subroutine csplinec(xi,yi,y,dydx1,dydx2)
    ! better way of doing things http://mathworld.wolfram.com/CubicSpline.html
    use lin_alg, only : tridiag
    !===============================================!
    ! Cubic spline                                  !
    ! Daniel Celis Garza 9 Apr. 2016                !
    !-----------------------------------------------!
    ! Inputs:                                       !
    ! dydx = derivatives                            !
    ! xi   = points in x                            !
    ! yi   = points in y                            !
    ! dxdy = array of derivatives at points 1 and N !
    !-----------------------------------------------!
    ! Locals:                                       !
    ! a() =
    ! b() =
    ! c() =                        !
    ! d() =
    !-----------------------------------------------!
    ! Outputs:                                      !
    ! y() = array of interpolated values            !
    !===============================================!
    implicit none
    real(dp), intent(in)           :: xi(:), yi(:)
    real(dp), optional, intent(in) :: dydx1, dydx2
    real(dp), intent(out)          :: y(:)
    real(dp), dimension(size(xi))  :: a, b, c, d
    integer                        :: n

    check_size: if (size(xi) /= size(yi) .or. size(xi) /= size(y)) then
      print*, " Error: subr: csplinec: check_size: sizes of xi, yi and y arrays must be the same. "
      return
    end if check_size

    ! Number of points.
    n = size(xi)

    ! Set up tridiagonal matrix.

    ! Upper diagonal.
    ! a(1) = x(2) - x(1)
    ! ...
    ! a(n-1) = x(n) - x(n-1)
    a(1: n - 1) = xi(2: n) - xi(1: n - 1)

    ! Right hand side (image of the matrix)
    ! b(1) = 6 * ( y(2) - y(1) ) / ( x(2) - x(1) )
    ! ...
    ! b(n-1) = 6 * ( y(n) - y(n-1) ) / ( x(n) - x(n - 1) )
    b(1: n - 1) = 6._dp * ( yi(2: n) - yi(1: n - 1) ) / a(1: n - 1)
    ! b(2) = 6 *( ( y(3) - y(2) ) / ( x(3) - x(2) ) - ( y(2) - y(1) ) / ( x(2) - x(1) ) )
    ! ...
    ! b(n-1) = 6 *( ( y(n) - y(n - 1) ) / ( x(n) - x(n - 1) ) - ( y(n - 1) - y(n - 2) ) / ( x(n - 1) - x(n - 2) ) )
    b(2: n - 1) = b(2: n - 1) - b(1: n - 2)

    ! Lower diagonal.
    ! c(2) = x(2) - x(1)
    ! ...
    ! c(n-1) = x(n - 1) - x(n - 2)
    c(2: n - 1) = a(1: n - 2)

    ! Main diagonal.
    ! d(2) = 2 * ( x(3) - x(2) + x(2) - x(1) )
    ! d(2) = 2 * ( x(3) - x(1) )
    ! ...
    ! d(n-1) = 2 * ( x(n - 1) - x(n - 2) + x(n - 2) - x(n - 3) )
    d(2: n - 1) = 2._dp * ( a(2: n - 1) + c(2: n - 1) )
    ! Because A = 1 at x1. numerical recipes 107--108
    d(1) = 1._dp
    ! Because B = 1 at xn. numerical recipes 107--108
    d(n) = 1._dp

    ! Lower boundary conditions.
    ! If the Lower Boundary Condition (LBC) is not present.
    lbc: if (present(dydx1) .eqv. .false.) then
      ! If there is no LBC, it is set to the natural value.
      a(1) = 0._dp
      b(1) = 0._dp
    else lbc
      ! Else, it is set to a specified first derivative.
      a(1) = 0.5_dp
      b(1) = ( 3._dp / ( xi(2) - xi(1) ) ) * ( ( yi(2) - yi(1) ) / ( xi(2) - xi(1) ) - dydx1 )
    end if lbc

    ! Upper boundary conditions.
    ! If the Upper Boundary Condition (UBC) is not present.
    ubc: if (present(dydx2) .eqv. .false.) then
      ! If there is no LBC, it is set to the natural value.
      c(n) = 0._dp
      b(n) = 0._dp
    else ubc
      ! Else, it is set to a specified first derivative.
      c(n) = 0.5_dp
      b(n) = ( - 3._dp / (http://www.paulinternet.nl/?page=bicubic xi(n) - xi(n - 1) ) ) * ( ( yi(n) - yi(n - 1) ) / ( xi(n) - xi(n - 1) ) - dydx2 )
    end if ubc

    ! Call subroutine to solve tridiagonal matrix.
    call tridiag(c(2: n), d, a(1: n - 1), b, y)
  end subroutine csplinec

  subroutine csplnei(xi, yi, yi2, x, y)
    implicit none
    real(dp), intent(in)  :: x, xi(:), yi(:), yi2(:)
    real(dp), intent(out) :: y
    real(dp) :: h, a, b
    integer  :: n, oidx, nidx

    ! Set the value of n to be the number of dots to interpolate.
    n = size(xi)
    ! Check that the sizes of all arrays involved are the same.
    check_size: if ( n /= size(yi) .or. n /= size(yi2) ) then
      write(*,*) " Error: subr: csplinei: check_size: sizes of xi, yi and yi2 arrays must be the same; size(xi) = ", n, &
      "; size(yi) = ", size(yi), "; size(yi2) = ", size(yi2)
      return
    end if check_size

    ! Find the index corresponding to the value of x.
    oidx = vbisect(xi, x)

    nidx = oidx + 1

    ! Calculate the interval.
    h = xi( nidx ) - xi( oidx )
    ! Check if h is zero.
    if ( h == 0._dp ) write(*,*) " Error: subr: csplinei: h is too small, wrong values of xi; h = ", h, "; x(nidx) = ", xi(nidx), &
    "; x(oidx) = ", xi(oidx)

    a = ( x - xi( oidx ) ) / h
    b = ( xi( nidx ) - x  ) / h

    ! Construct the cubic equation.
    y = a * yi(nidx) + b * yi(oidx) + ( ( a*a*a - a ) * yi2(nidx) + ( b*b*b - b ) * yi2(oidx) ) * h*h / 6._dp
  end subroutine csplnei

  function vbisect(xi, x)
    ! Find x in within an interval in xa.
    implicit none
    real(dp), intent(in) :: xi(:), x
    integer :: vbisect
    integer :: n, i, a, b, c
    logical :: up

    n  = size(xi)
    up = ( xi(n) >= xi(1) )
    a  = 0
    b  = n + 1

    ! Carry out Bisection.
    bisect: do
      ! If we are in adjacent or equal positions, exit.
      ! [ | n | n+1 | | ]
      ! Exit if (a = n and b = n + 1), or  Exit if ( a = n and b = n )
      if ( b - a <= 1 ) exit bisect
      ! Find midpoint via integer division.
      c = ( a + b ) / 2
      ! Decide whether to move a Up or b Down.
      ud: if ( up .eqv. ( x >= xi(c) ) ) then
        a = c
      else ud
        b = c
      end if ud
    end do bisect

    ! Check Location of X Within the Interval (LXWI).
    lxwi: if ( x == xi(1) ) then
      vbisect = 1
    else if ( x == xi(n) ) then
      vbisect = n - 1
    else lxwi
      vbisect = a
    end if lxwi
  end function vbisect

  function cctrilinint(grid, x, idx, grid_vol)
    ! http://bmia.bmt.tue.nl/people/BRomeny/Courses/8C080/Interpolation.pdf
    ! There is a typo in the pdf. In trilinear interpolation formula z_0 -> z_8 should be v_0 -> v8.
    ! Cubic Cell TRILINnear INTerpolation.
    ! Uses a cubic cell to carry out trilinear interpolation.
    implicit none
    real(dp), intent(in) :: grid(:,:,:), x(:), grid_vol ! 3D grid, value of x, grid volume.
    integer, intent(in)  :: idx(:) ! Grid indices.
    real(dp) :: cctrilinint ! Value of the trilinear interpolation.
    integer :: idxp1(size(idx))

    idxp1 = idx + 1
                  ! grid(0,0,0)
    cctrilinint = grid(  idx(1),   idx(2),   idx(3)) * product(idxp1 - x) + &
                  ! grid(1,0,0)
                  grid(idxp1(1),   idx(2),   idx(3)) * (- idx(1) + x(1)) * (idxp1(2) - x(2)) * (idxp1(3) - x(3)) + &
                  ! grid(0,1,0)
                  grid(  idx(1), idxp1(2),   idx(3)) * (idxp1(1) - x(1)) * (- idx(2) + x(2)) * (idxp1(3) - x(3)) + &
                  ! grid(1,1,0)
                  grid(idxp1(1), idxp1(2),   idx(3)) * (- idx(1) + x(1)) * (- idx(2) + x(2)) * (idxp1(3) - x(3)) + &
                  ! grid(0,0,1)
                  grid(  idx(1),   idx(2), idxp1(3)) * (idxp1(1) - x(1)) * (idxp1(2) - x(2)) * (- idx(3) + x(3)) + &
                  ! grid(1,0,1)
                  grid(idxp1(1),   idx(2), idxp1(3)) * (- idx(1) + x(1)) * (idxp1(2) - x(2)) * (- idx(3) + x(3)) + &
                  ! grid(0,1,1)
                  grid(  idx(1), idxp1(2), idxp1(3)) * (idxp1(1) - x(1)) * (- idx(2) + x(2)) * (- idx(3) + x(3)) + &
                  ! grid(1,1,1)
                  grid(idxp1(1), idxp1(2), idxp1(3)) * product(- idx + x)

    cctrilinint = cctrilinint/grid_vol
  end function cctrilinint

end module interpolation
