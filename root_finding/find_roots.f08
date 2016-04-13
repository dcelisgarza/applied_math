module find_roots
  use nrtype
contains
  subroutine outbracket(func,interval,success,max_it,factor)
    implicit none
    interface cont_func
      function func(x)
        use nrtype
        implicit none
        real(dp), intent(in)  :: x
        real(dp)              :: func
      end function func
    end interface cont_func
    real(dp), intent(inout)   :: interval(2)
    real(dp), intent(in)      :: factor
    integer, intent(in)       :: max_it
    logical, intent(out)      :: success
    integer                   :: i
    real(dp)                  :: f1, f2

    failsafedo: do
      failsafeif: if (interval(1) == interval(2)) then
        write(*,*) 'Please provide a proper interval [A, B] where A < B.'
        write(*,'(A)', advance = 'no') ' A = '
        read*, interval(1)
        write(*,'(A)', advance = 'no') ' B = '
        read*, interval(2)
      else failsafeif
        exit failsafedo
      end if failsafeif
    end do failsafedo

    f1 = func(interval(1))
    f2 = func(interval(2))
    success = .true.
    expand_int: do i = 1, max_it
      if ( f1 * f2 < 0.) return
      expand: if ( abs(f1) < abs(f2) ) then
        interval(1)  = interval(1) + factor * (interval(1) - interval(2))
        f1 = func(interval(1))
      else expand
        interval(2)  = interval(2) + factor * (interval(2) - interval(1))
        f2 = func(interval(2))
      end if expand
    end do expand_int
    success = .false.
  end subroutine outbracket

  subroutine inbracket(func,mast_interval,interval_seg,n_seg,n_root)
    implicit none
    interface cont_func
      function func(x)
        use nrtype
        implicit none
        real(dp), intent(in)  :: x
        real(dp)              :: func
      end function func
    end interface cont_func
    real(dp), intent(in)    :: mast_interval(2)
    integer, intent(in)     :: n_seg
    integer, intent(inout)  :: n_root ! # of segments, # of roots
    real(dp), intent(out)   :: interval_seg(:,:)
    real(dp)                :: fa, fb, x, dx
    integer                 :: i, int_count

    int_count = 0
    x  = mast_interval(1)
    dx = (mast_interval(2) - mast_interval(1)) / n_seg
    fa = func(x)
    interval_loop: do i = 1, n_seg
      ! Traveling right along the x axis, because we want to calculate B of an internal interval.
      x  = x + dx
      fb = func(x)
      check_sign: if (fa * fb < 0.) then
        int_count = int_count + 1
        interval_seg(1,int_count) = x - dx
        interval_seg(2,int_count) = x
        if (int_count == n_root) exit interval_loop
      end if check_sign
      fa = fb
    end do interval_loop
    if (i-1 == n_seg) write(*,*) ' Overflowed number of segments, add more.'
    n_root = int_count
  end subroutine inbracket

  subroutine bisection(func,interval,tol,max_it,root)
    ! Bisection subroutine.
    implicit none
    interface cont_func
      function func(x)
        use nrtype
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: func
      end function func
    end interface cont_func
    real(dp), intent(in)     :: interval(2), tol ! Interval, tolerance
    integer, intent(in)      :: max_it ! Maximum # of iterations
    real(dp), intent(out)    :: root
    integer                  :: i ! Counter.
    real(dp)                 :: a, b, c, fa, fc

    a = interval(1)
    b = interval(2)
    fa = func(a)
    fc = func(c)

    bisecting: do i = 1, max_it
      c  = .5_dp*(a+b)
      fc = func(c)
      if ( fc == 0._dp .or. .5_dp*(b - a) < tol) exit bisecting
      check_sign: if ( fc * fa < 0._dp ) then
        b = c
      else check_sign
        a  = c
        fa = func(a)
      end if check_sign
    end do bisecting

    found_root: if (i > max_it) then
      write(*,*) ' Iterations not sufficient for root finding.'
    else found_root
      root = c
    end if found_root
  end subroutine bisection
end module find_roots
