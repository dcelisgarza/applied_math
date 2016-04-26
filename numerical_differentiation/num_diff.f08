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
      dncoef1(n+1-i) = (-1)**i * rcsv_bicoef(n,i)
    end do cc
  end function dncoef1

  function dncoefn(m, n, a, xi) result(rdncoefn)
    ! n = order of the derivative.
    ! k = number of grid points (h-steps)
    ! xi
    ! rdncoefn = array of coefficients of nth order with k - n + 1 order of accuracy
    ! http://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0935077-0/S0025-5718-1988-0935077-0.pdf
    implicit none
    integer, intent(in)            :: m, n
    real(dp), intent(in)           :: xi, a(0:n)
    !real(dp), optional             :: xi
    real(dp)                       :: rdncoefn(0:n)
    integer  :: sn, nu, sm
    real(dp) :: c1, c2, c3
    real(dp) :: d(0:m,0:n,0:n)

    ! Check whether we give xi, if not default to xi = 0
    !if(.not. present(xi)) xi = 0._dp

    ! Check array size.
    !check_size: if( size(a) /= size(rdncoefn) .or. size(a) /= m + 1 ) then
    !  write(*,*) " Error: num_diff: dncoefn: check_size: Grid size and the size of the coefficient table are not equal to n + 1. &
    !               size(grid) = ", size(a), " size(rdncoefn) = ", size(rdncoefn), " n + 1 = ", n + 1
    ! return
    !end if check_size
    d = 0._dp
    d(0,0,0)  = 1._dp
    c1 = 1._dp

    ! Loop Over Derivatives.
    od: do sn = 1, n

      c2 = 1._dp

      ! Loop Over H-Steps.
      ohs: do nu = 0, sn - 1
        c3 = a(sn) - a(nu)
        c2 = c2 * c3

        ! Inner loop Over Derivatives.
        iod: do sm = 0, minval([sn, m])
          ! Prevent Memory Dump
          pmd1: if( sm /= 0 ) then
            d(sm, sn, nu) = ( a(sn) - xi ) * d(sm,sn-1,nu) - sm * d(sm-1,sn-1,nu)
          else pmd1
            d(sm, sn, nu) = ( a(sn) - xi ) * d(sm,sn-1,nu)
          end if pmd1
          d(sm, sn, nu) = d(sm, sn, nu) / c3
        end do iod
      end do ohs

      ! loop over derivatives
      do sm = 0, minval([sn,m])
        pmd2: if ( sm /= 0 ) then
          d(sm, sn, sn) = c1/c2 * ( sm * d(sm-1,sn-1,sn-1) - (a(sn-1) - xi) *   d(sm,sn-1,sn-1) )
        else pmd2
          d(sm, sn, sn) = -c1/c2 * (a(sn-1) - xi) *   d(sm,sn-1,sn-1)
        end if pmd2
      end do
      c1 = c2

    end do od
    print*, d(2,4,:)
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
      fdo1 = fdo1 + coef1(i+1) * func( x - (n-i)*h )
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
      bdo1 = bdo1 + coef1(n-i+1) * func( x - i*h )
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
        cdo1 = cdo1 + coef1(n-i+1) * func( x + (n/2-i)*h )
      end do cde
    else cne
      cdo: do i = 0, n
        cdo1 = cdo1 + coef1(n-i+1) * ( func( x - h/2._dp + (n/2._dp-i)*h ) + func( x + h/2._dp + (n/2._dp-i)*h ) )/2_dp
      end do cdo
    end if cne
    cdo1 = cdo1/h**n
  end function cdo1

end module num_diff
