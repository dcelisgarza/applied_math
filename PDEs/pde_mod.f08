module finite_elements
  use nrtype
contains

  subroutine boundry_val(derivs,x,yi,yo,num_diff,h)
    implicit none

    real(dp), intent(in)  :: x(:), yi(:), h(:)
    real(dp), intent(out) :: yo(:)

    interface numerical_differentiation
      function num_diff(x,y,h)
        use nrtype
        implicit none
        real(dp), intent(in) :: x, h
        real(dp), intent(in) :: y
        real(dp)             :: num_diff(:)
      end function num_diff
    end interface numerical_differentiation

    !interface derivatives
    !    subroutine derivs(x,y,dydx)
    !      use nrtype
    !      implicit none
    !      real(dp), intent(in)  :: x(:), y(:)
    !      real(dp), intent(out) :: dydx(:)
    !  end subroutine derivs
    !end interface derivatives
    integer               :: i


  end subroutine boundry_val

end module finite_elements
