module num_diff
  use nrtype
  use nrtools, only : rcsv_bicoef
contains
  function sifodi1(func,x,h)
    ! si := simple, fo := forward, di1 = 1st derivative
    implicit none
    real(dp), intent(in) :: x, h ! Variable and step size.
    real(dp)             :: sifodi1
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    div0: if (h == 0) then
      write(*,*) ' Error: sifodi1:: Division by zero.'
      return
    end if div0
    sifodi1 = (func(x+h) - func(x))/h
  end function sifodi1

  function eicedi1(func,x,h)
    ! ei := eigth order, ce := central, di1 := 1st derivative
    implicit none
    real(dp), intent(in) :: x, h ! Variable and step size.
    real(dp)             :: eicedi1
    real(dp)             :: ho8, ho8p6   ! h/8
    real(dp), parameter  :: c4 = 1./280, c3 = 4./105, c2 = 0.2, c1 = 0.8
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    ho8   = 0.125_dp*h
    ho8p6 = ho8**6._dp

    div0: if (ho8p6 == 0) then
      write(*,*) ' Error: eicedi1 :: Division by zero.'
      return
    end if div0
    eicedi1 = (c4*func(x-ho8*4._dp)-c3*func(x-ho8*3._dp)+c2*func(x-ho8*2._dp)&
    -c1*func(x-ho8)&
    +c1*func(x+ho8)&
    -c4*func(x+ho8*4._dp)+c3*func(x+ho8*3._dp)-c2*func(x+ho8*2._dp))/ho8p6
  end function eicedi1

  function dncoef1(n)
    ! Coefficients for n'th order derivatives with 1st order errors for forward and backward differences, 2nd order errors for central differences
    ! https://en.wikipedia.org/wiki/Finite_difference
    implicit none
    integer, intent(in) :: n
    real(dp)            :: dncoef1(n+1)
    integer :: i

    ! Calculate the Coefficients.
    cc: do i = 0, n
      dncoef1(i+1) = (-1)**i * float(rcsv_bicoef(n,i))
    end do cc
  end function dncoef1

  function dncoefn(n,k)
    ! TO DO:
    ! Make this a recursive function or a loop.


    ! Coefficients for n'th order derivatives with k'th order errors for forward and backward, 2k'th order errors for central differences
    ! https://en.wikipedia.org/wiki/Finite_difference
    ! https://en.wikipedia.org/wiki/Finite_difference_coefficient
    ! https://en.wikipedia.org/wiki/Taylor_series
    implicit none
    integer, intent(in) :: n, k
    real(dp)            :: dncoefn(n+k), dkcoefk(n+k+2)
    integer :: c

    if ( k < n ) write(*,*) " Error: dncoefn: k must be greater than n, k = ", k, " n = ", n

    ! Calculate the difference between the order of the error and order of the derivative.
    c = k - n

    ! Calculate the coefficients for the n'th derivative.
    dkcoefk(1:n+1) = dncoef1(n)
    !print*, dkcoefk(1:n+1)

    ! Calculate the coefficients for the k'th order error.
    dkcoefk(n+2:n+k+2) = - dncoef1(k) / k

    ! Map coefficients.
    dncoefn(1:c)     = dkcoefk(n+2:n+1+c)
    dncoefn(c+1:n+k) = dkcoefk(1:n+1) + dkcoefk(n+2+c:n+k+2)

  end function dncoefn

  function fdo1(func, x, h, coef1)
    ! forward numerical differentiator for n'th derivative with first order errors
    implicit none
    real(dp), intent(in) :: x, h
    integer, intent(in)  :: coef1(:)
    real(dp)             :: fdo1
    integer :: i, n
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    fdo1 = 0._dp
    ! Size of the coefficient array.
    n = size(coef1)-1
    ! Calculate Derivative with .
    cd: do i = 0, n
      fdo1 = fdo1 + coef1(n-i+1) * func( x - (n-i)*h )
    end do cd
    fdo1 = fdo1/h**n
  end function fdo1

  function bdo1(func, x, h, coef1)
    ! backward numerical differentiator for n'th derivative with first order errors
    implicit none
    real(dp), intent(in) :: x, h
    integer, intent(in)  :: coef1(:)
    real(dp)             :: bdo1
    integer :: i, n
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    bdo1 = 0._dp
    ! Size of the coefficient array.
    n = size(coef1)-1
    ! Calculate Derivative with .
    cd: do i = 0, n
      bdo1 = bdo1 + coef1(i+1) * func( x - i*h )
    end do cd
    bdo1 = bdo1/h**n
  end function bdo1

  function cdo1(func, x, h, coef1)
    ! central numerical differentiator for n'th derivative with first order errors
    implicit none
    real(dp), intent(in) :: x, h
    integer, intent(in)  :: coef1(:)
    real(dp)             :: cdo1
    integer :: i, n
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    cdo1 = 0._dp
    ! Size of the coefficient array.
    n = size(coef1)-1
    ! Check if N is Even.
    cne: if (mod(n,2) == 0) then
      cde: do i = 0, n
        cdo1 = cdo1 + coef1(i+1) * func( x + (n/2-i)*h )
      end do cde
    else cne
      cdo: do i = 0, n
        cdo1 = cdo1 + coef1(i+1) * ( func( x - h/2._dp + (n/2._dp-i)*h ) + func( x + h/2._dp + (n/2._dp-i)*h ) )/2_dp
      end do cdo
    end if cne
    cdo1 = cdo1/h**n
  end function cdo1

end module num_diff
