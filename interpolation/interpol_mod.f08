module interpolation
  use nrtype
contains
  subroutine poly_interpol(xi,yi,x,y,dy)
    implicit none
    real(dp), intent(in)          :: x, xi(:), yi(:)
    real(dp), intent(out)         :: y, dy
    integer                       :: i, j, idx, aux_idx(1), n
    real(dp), dimension(size(xi)) :: c, d, difx, aux_difx, den
    ! Making sure we're using equally sized arrays.
    if(size(xi) /= size(yi)) then
      print*, ' Size of x and y arrays must be equal.'
      return
    endif
    ! The order of the interpolation.
    n = size(xi)
    ! Initialising relationships of Cs and Ds.
    c = yi
    d = yi
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

    ! Looping over the columns of Neville's tree.
    nev_loop: do i = 1, n - 1
      den(1:n-i) = difx(1:n-i) - difx(1+i:n)
      if ( any( den(1:n-i) == 0.0_dp ) ) print*, ' Error: subr:: poly_interpol: nev_loop. Division by zero.'
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

  end subroutine poly_interpol
end module interpolation
