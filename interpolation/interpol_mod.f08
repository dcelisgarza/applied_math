module interpolation
  use numbers
contains
  subroutine poly_interpol(xi,yi,n,x,y)
    implicit none
    real(dp), intent(in)  :: xi(:), yi(:)
    integer, intent(in)   :: n
    real(dp), intent(out) :: x, y
    real(dp)              :: dif1, dif2
    integer               :: i, idx
    real(dp)              :: c(size(yi)), d(size(yi))



    ! Some sort of recursive relationship
    idx  = 1

    dif1 = abs(x-xi(1))
    lclosest_index:do i = 1, n
      dif2 = abs(x-xi(i))
      iclosest_index: if (dif2 < dif1) then
        idx  = i
        dif1 = dif2
      end if iclosest_index
      c(i) = yi(i)
      d(i) = yi(i)
    end do lclosest_index
    y   = y(idx)
    idx = idx - 1

  end subroutine poly_interpol
end module interpolation
