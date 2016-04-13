module interpolation
  use lin_alg
contains
  subroutine pinex(xi,yi,x,y,dy)
    ! Polynomial INterpolation and EXtrapolation (PINEX).
    implicit none
    real(dp), intent(in)          :: x, xi(:), yi(:)
    real(dp), intent(out)         :: y, dy
    integer                       :: i, j, idx, aux_idx(1), n
    real(dp), dimension(size(xi)) :: c, d, difx, aux_difx, den

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

      err: if (2*idx < n-i) then
        dy  = c(idx+1)
      else err
        dy  = d(idx)
        idx = idx-1
      end if err
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
        print*, " Size of x and y arrays must be equal."
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

  subroutine cspline(xi,yi,y,dydx1,dydx2)
    !===============================================!
    ! Cubic spline                                  !
    ! Runge-Kutta 4 Gill's method                   !
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
      print*, " Error: subr: cspline: check_size: sizes of xi, yi and y arrays must be the same."
      return
    end if check_size

    ! Number of points.
    n = size(xi)

    ! Set up tridiagonal matrix.

    ! a(1) = x(2) - x(1)
    ! ...
    ! a(n-1) = x(n) - x(n-1)
    a(1: n - 1) = xi(2: n) - xi(1: n - 1)

    ! b(1) = 6 * ( y(2) - y(1) ) / ( x(2) - x(1) )
    ! ...
    ! b(n-1) = 6 * ( y(n) - y(n-1) ) / ( x(n) - x(n - 1) )
    b(1: n - 1) = 6._dp * ( yi(2: n) - yi(1: n - 1) ) / a(1: n - 1)

    ! b(2) = 6 *( ( y(3) - y(2) ) / ( x(3) - x(2) ) - ( y(2) - y(1) ) / ( x(2) - x(1) ) )
    ! ...
    ! b(n-1) = 6 *( ( y(n) - y(n - 1) ) / ( x(n) - x(n - 1) ) - ( y(n - 1) - y(n - 2) ) / ( x(n - 1) - x(n - 2) ) )
    b(2: n - 1) = b(2: n - 1) - b(1: n - 2)

    ! c(2) = x(2) - x(1)
    ! ...
    ! c(n-1) = x(n - 1) - x(n - 2)
    c(2: n - 1) = a(1: n - 2)

    ! d(2) = 2 * ( x(3) - 2 * x(2) + x(1) )
    ! ...
    ! d(n-1) = 2 * ( x(n) - 2 * x(n - 1) + x(n - 2) )
    d(2: n - 1) = 2._dp * ( a(2: n - 1) - c(2: n - 1) )
    d(1) = 1._dp
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
      b(1) = ( 3._dp / ( xi(2) - x(1) ) ) * ( ( yi(2) - yi(1) ) / ( xi(2) - xi(1) ) - dydx1 )
    end if lbc

    ! Lower boundary conditions.
    ! If the Lower Boundary Condition (LBC) is not present.
    lbc: if (present(dydx1) .eqv. .false.) then
      ! If there is no LBC, it is set to the natural value.
      a(1) = 0._dp
      b(1) = 0._dp
    else lbc
      ! Else, it is set to a specified first derivative.
      c(1) = 0.5_dp
      b(1) = ( - 3._dp / ( xi(n) - x(n - 1) ) ) * ( ( yi(n) - yi(n - 1) ) / ( xi(n) - xi(n - 1) ) - dydx2 )
    end if lbc

    ! Call subroutine to solve tridiagonal matrix.
     call tridiag(d(2: n), c, a(1: n - 1), b, y)
  end subroutine cspline
end module interpolation
