module find_roots
  use numbers
contains
  subroutine bisection(func,interval,tol,max_it,root)
    implicit none
    interface cont_func
      function func(x)
        use numbers
        implicit none
        real(dp), intent(in),optional :: x
        real(dp)             :: func
      end function func
    end interface cont_func
    real(dp), intent(in) :: interval(2), tol ! Interval, tolerance
    integer, intent(in)  :: max_it ! Maximum # of iterations
    real(dp), intent(out):: root
    integer              :: i ! Counter.
    real(dp)             :: a, b, c

    a = interval(1)
    b = interval(2)

    bisecting: do i = 1, max_it
      c = .5*(a+b)
      if ( func(c) == 0._dp .or. .5_dp*(b - a) < tol) exit bisecting
      check_sign: if ( (func(c) < 0._dp .and. func(a) > 0._dp) .or. (func(c) > 0._dp .and. func(a) < 0._dp) ) then
        b = c
      else check_sign
        a = c
      end if check_sign
    end do bisecting
    root = c
    end subroutine bisection
end module find_roots
