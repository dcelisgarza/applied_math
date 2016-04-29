module find_roots
  !=========================================================!
  ! Module for root finding.                                !
  ! Author: Daniel Celis Garza                              !
  ! Last modified: April 28, 2016                           !
  ! Contains:                                               !
  ! subroutine outbracket: Expands the interval until it    !
  ! crosses a root (function changes sign).                 !
  ! subroutine inbracket: Surrounds roots within intervals  !
  ! subroutine bisection: Carries out bisection to find a   !
  ! root within an interval.                                !
  !=========================================================!
  use nrtype
contains
  subroutine outbracket(func, interval, max_it, success, factor)
    !=======================================================!
    ! Gradually expands bracket until it finds a root.      !
    !-------------------------------------------------------!
    ! Cycles through a number of iterations until the       !
    ! maximum number of iterations is reached. The interval !
    ! is expanded according to a user-provided factor. If   !
    ! the cycle finds a root, it exits the subroutine with  !
    ! success = .true., else it exits with success = .false.!
    ! and a greater interval.                               !
    !-------------------------------------------------------!
    ! Inputs:                                               !
    ! func   = function whose roots we're trying to bracket !
    ! factor = optional factor for increasing the interval  !
    ! max_it = maximum number of iterations                 !
    !-------------------------------------------------------!
    ! Outputs:                                              !
    ! success = logical value stating whether a root was    !
    !           bracketted                                  !
    !-------------------------------------------------------!
    ! Inputs-Outputs:                                       !
    ! interval() = input a proposed interval, output an     !
    !              interval which brackets the root.        !
    !-------------------------------------------------------!
    ! Locals:                                               !
    ! f() = value of func() at the interval bounds          !
    ! i   = counter                                         !
    !=======================================================!
    implicit none
    interface cont_func
      function func(x)
        use nrtype
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: func
      end function func
    end interface cont_func
    real(dp), intent(inout) :: interval(2)
    integer,  intent(in)    :: max_it
    logical,  intent(out)   :: success
    real(dp), optional      :: factor
    real(dp)                :: f(2)
    integer                 :: i

    ! Check whether factor is given, if it's not give it a default value of 2.
    if( .not. present(factor) ) factor = 2._dp

    ! Loop to check whether the Interval given is Appropriate and ask the user to provide an appropriate one.
    iado: do
      iaif: if ( interval(1) == interval(2) ) then
        write(*,*) " Please provide a proper interval [A, B] where A < B."
        write(*,"(A)", advance = "no") " A = "
        read*, interval(1)
        write(*,"(A)", advance = "no") " B = "
        read*, interval(2)
      else iaif
        exit iado
      end if iaif
    end do iado

    ! Evaluate the lower and upper bounds of the interval on the fuction provided.
    f(1) = func( interval(1) )
    f(2) = func( interval(2) )

    ! Define success as true in case the subroutine finds a root.
    success = .true.

    ! Loop to EXPand the INTerval
    expintdo: do i = 1, max_it
      if ( f(1) * f(2) < 0.) return
      ! EXPand the INTerval
      expintif: if ( abs(f(1)) < abs(f(2)) ) then
        interval(1)  = interval(1) + factor * (interval(1) - interval(2))
        f(1)         = func(interval(1))
      else expintif
        interval(2)  = interval(2) + factor * (interval(2) - interval(1))
        f(2)         = func(interval(2))
      end if expintif
    end do expintdo

    ! If the subroutine fails to find a root, succes is false.
    success = .false.

  end subroutine outbracket

  subroutine inbracket(func, intvl, n_seg, sub_intvl, n_root)
    !=======================================================!
    ! Find intervals which enclose a given # of roots       !
    !-------------------------------------------------------!
    ! Cycles through a number of iterations until the       !
    ! number of proposed roots is reached. A master interval!
    ! is traversed, and subintervals are generated when a   !
    ! root is found (function changes sign). Outputs the    !
    ! intervals found and the number of roots found.        !
    !-------------------------------------------------------!
    ! Inputs:                                               !
    ! func  = function whose roots we're trying to bracket  !
    ! intvl = master interval (contains the subints)        !
    ! n_seg = number of segments to analyse                 !
    !-------------------------------------------------------!
    ! Outputs:                                              !
    ! sub_intvl = subintervals                              !
    !-------------------------------------------------------!
    ! Inputs-Outputs:                                       !
    ! n_root() = # of roots found                           !
    !-------------------------------------------------------!
    ! Locals:                                               !
    ! f()     = value of func() at the subinterval bounds   !
    ! x       = independent variable                        !
    ! dx      = step size                                   !
    ! i       = counter                                     !
    ! n_intvl = # of intervals                              !
    !=======================================================!
    implicit none
    interface cont_func
      function func(x)
        use nrtype
        implicit none
        real(dp), intent(in) :: x
        real(dp)             :: func
      end function func
    end interface cont_func
    real(dp), intent(in)    :: intvl(2)
    integer,  intent(in)    :: n_seg
    real(dp), intent(out)   :: sub_intvl(:,:)
    integer,  intent(inout) :: n_root
    real(dp)                :: f(2), x, dx
    integer                 :: i, n_intvl

    ! Initialise the number of intervals found.
    int_count = 0

    ! Set initial x to the value of the master interval.
    x  = intvl(1)
    ! Calculate the step size from the number of desired segments and master interval.
    dx = (intvl(2) - intvl(1)) / n_seg

    ! Evaluate function at lower end of the master interval.
    f(1) = func(x)

    ! Find Root Bracketting Intervals.
    frbido: do i = 1, n_seg
      ! Travel right along the x axis, because we want to calculate B of an internal interval.
      x  = x + dx
      ! Evaluate function at prospective upper bound for the interval.
      f(2) = func(x)
      ! Check whether we have found a root.
      frbiif: if ( f(1) * f(2) < 0. ) then
        ! If a root has been found, increase the interval counter.
        n_intvl = n_intvl + 1
        ! Record the lower and upper bounds of the subinterval.
        sub_intvl(1,n_intvl) = x - dx
        sub_intvl(2,n_intvl) = x
        ! If we have reached the number of roots that we want to bracket, exit subroutine.
        if (n_intvl == n_root) exit frbido
      end if frbiif
      ! Set up for the next subinterval.
      f(1) = f(2)
    end do frbido

    ! Check if the master interval is big enough to find all proposed roots.
    if ( i-1 == n_seg ) write(*,*) " Increase range of the master interval or change the number of segments."
    ! Print the roots found.
    n_root = n_intvl
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
      write(*,*) " Iterations not sufficient for root finding."
    else found_root
      root = c
    end if found_root
  end subroutine bisection
end module find_roots
