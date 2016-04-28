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

  function dncoefn(n, k, grid, xi) result(coefnk)
    ! n = order of the derivative.
    ! k = number of grid points (h-steps)
    ! xi
    ! rdncoefn = array of coefficients of nth order with k - n + 1 order of accuracy
    ! http://www.ams.org/journals/mcom/1988-51-184/S0025-5718-1988-0935077-0/S0025-5718-1988-0935077-0.pdf
    implicit none
    integer, intent(in)    :: n, k
    real(dp), intent(in)   :: xi, grid(0:k)
    integer  :: sk, nu, sn
    real(dp) :: c1, c2, c3
    real(dp) :: coefnk(0:n,0:k,0:k)

    ! Check for parameter errors.

    coefnk = 0._dp
    coefnk(0,0,0)  = 1._dp
    c1 = 1._dp

    ! Loop Over Derivatives.
    od: do sk = 1, k

      c2 = 1._dp

      ! Loop Over H-Steps.
      ohs: do nu = 0, sk - 1
        c3 = grid(sk) - grid(nu)
        c2 = c2 * c3

        ! Inner loop Over Derivatives.
        iod: do sn = 0, minval([sk, n])
          ! Prevent Memory Dump
          pmd1: if( sn /= 0 ) then
            coefnk(sn, sk, nu) = ( grid(sk) - xi ) * coefnk(sn,sk-1,nu) - sn * coefnk(sn-1,sk-1,nu)
          else pmd1
            coefnk(sn, sk, nu) = ( grid(sk) - xi ) * coefnk(sn,sk-1,nu)
          end if pmd1
          coefnk(sn, sk, nu) = coefnk(sn, sk, nu) / c3
        end do iod
      end do ohs

      ! loop over derivatives
      do sn = 0, minval([sk,n])
        pmd2: if ( sn /= 0 ) then
          coefnk(sn, sk, sk) = c1/c2 * ( sn * coefnk(sn-1,sk-1,sk-1) - (grid(sk-1) - xi) *   coefnk(sn,sk-1,sk-1) )
        else pmd2
          coefnk(sn, sk, sk) = -c1/c2 * (grid(sk-1) - xi) *   coefnk(sn,sk-1,sk-1)
        end if pmd2
      end do
      c1 = c2
    end do od
  end function dncoefn

  function ndifnk(func, x, h, n, grid, coefn)
    ! forward/backward numerical differentiator for n'th derivatives with k'th order errors
    implicit none
    real(dp), intent(in) :: x, h, grid(:), coefn(:)
    integer, intent(in)  :: n
    real(dp)             :: ndifnk
    integer :: i, m
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in) :: x    ! Variables
        real(dp)             :: func ! Function
      end function func
    end interface functn

    ! Size of grid and coefn
    ! Check Size.
    check_size: if ( size(grid) == size(coefn) ) then
      ! Asign m to the size of the array.
      m = size(grid)
    else check_size
      write(*,*) " Error: num_diff: fdon: check_size: size(grid) must be equal to size(coefn); &
                   size(grid) = ", size(grid), " size(coefn) = ", size(coefn)
      return
    end if check_size

    ! Set fdon to 0.
    ndifnk = 0._dp

    ! Construct Derivative.
    cd: do i = 1, m
      ndifnk = ndifnk + coefn(i) * func( x + h * grid(i) )/h**n
    end do cd
  end function ndifnk

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
    ! central numerical differentiator for n'th derivatives with first order errors
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
