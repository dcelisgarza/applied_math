module interpolation
  use nrtype
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
      if ( any( den(1:n-i) == 0._dp ) ) print*, ' Error: subr:: poly_interpol: nev_loop. Division by zero.'
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
          print*, " Error:: rinex:: nev_loop: Denominator = 0, use different points."
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
end module interpolation
