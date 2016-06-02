module nrtools
  use nrtype
contains
  recursive function rcsv_bicoef(n, k) result(bicoef)
    !===========================================================================!
    ! Recursive function for the calculation of binomial coefficients. !
    ! April 19th 2016. !
    ! Daniel Celis Garza !
    !---------------------------------------------------------------------------!
    ! This function recursively calculates the binomial coefficient for a given !
    ! n and k. !
    ! n choose k = n!/( k! * (n - k)! ) = n choose (n - k) !
    !
    ! n!/( k! * (n - k)! ) = n(n-1)(n-2)...(n-k+1)/k                            !
    ! = (n/k) * (n-1)!/( (k-1)! * (n-k)! )                              !
    !---------------------------------------------------------------------------!
    ! Inputs:
    ! n = integer n
    ! k = integer k
    !---------------------------------------------------------------------------!
    ! Outputs:
    ! bicoef = integer type i16 !
    !===========================================================================!
    implicit none
    integer, intent(in) :: n, k
    real(dp)            :: bicoef

    ! Identify Key Scenarios.
    iks: if (k > n) then
      write(*,*) " Error: nrtools: rcsv_bicoef: iks: n must be greater than k. k = ", k, ", n = ", n
      return
    else if (k == 0._dp) then iks
      bicoef = 1._dp
    else if (k > n/2) then iks
      ! Recurse function so as to reduce the number of operations to at most n/2 at the numerator and denominator.
      bicoef = rcsv_bicoef(n, n-k)
    else iks
      ! Calculate binomial coefficient.
      bicoef = float(n)/float(k) * rcsv_bicoef(n-1, k-1)
    end if iks

  end function rcsv_bicoef

  function ln_gamma(x)
    ! Approximates the natural logarithm of the gamma function for x \in \reals and x > 0
    implicit none
    real(dp), intent(in) :: x
    real(dp) :: ln_gamma
    real(dp) :: tmp
    real(dp), parameter, dimension(7) :: coef = [1.000000000190015_dp,76.18009172947146_dp,-86.50532032941677_dp,&
    0._dp,0._dp,0._dp,0._dp]
    real(dp), parameter :: sqrtau = sqrt(tau)

    ! Check that x > 0.
    xg0: if (x <= 0._dp) then
      write(*,*) " Error: nrtools: ln_gamma: xg0: x must be greater than 0. x = ", x
      return
    end if xg0

    tmp = x + 5.5_dp
    tmp = (x + 0.5_dp) * log(tmp) - tmp
    ln_gamma = tmp + log(sqrtau * (coef(1) + sum( coef(2:7) / arthm_prog( x + 1._dp, 1._dp, size(coef(2:7)) ) ) ) / x )

  end function ln_gamma

  function arthm_prog(first, step, n)
    ! Calculates the arithmetic progression with step = step.
    implicit none
    real(dp), intent(in) :: first, step
    integer, intent(in)  :: n
    real(dp) :: arthm_prog(n)
    integer  :: i ! Counter

    arthm_prog(1) = first

    do i = 2, n
        arthm_prog(i) = arthm_prog(i-1) + step
    end do

  end function arthm_prog

  recursive function rcsv_gamma(z) result(gamma)

    ! Approximates the gamma function via Lanczos' approximation, with gamma = 8, N = 9
    ! https://en.wikipedia.org/wiki/Lanczos_approximation
    implicit none
    complex(dpc), intent(in) :: z
    complex(dpc) :: gamma
    complex(dpc) :: tmp, t
    real(dp) :: x
    integer  :: i
    real(dp), parameter, dimension(9) :: coef = [0.99999999999980993_dp, 676.5203681218851_dp, -1259.1392167224028_dp, &
    771.32342877765313_dp, -176.61502916214059_dp, 12.507343278686905_dp,&
    -0.13857109526572012_dp, 9.9843695780195716d-6, 1.5056327351493116d-7]
    real(dp), parameter :: sqrtau = sqrt(tau)

    tmp = z

    ! Calculate Gamma Function
    cgf: if( real(tmp) < 0.5_dp ) then

      gamma = pi / ( sin(pi*tmp) * rcsv_gamma(1._dp-tmp) )

    else cgf

      tmp = tmp - 1._dp
      x   = coef(1)

      ! Calculate the Arithmetic Progression and Divide.
      apd: do i = 2, size(coef)
        x = x + coef(i) / (tmp + i - 1)
      end do apd

      t = tmp + size(coef) - 1.5_dp

      gamma = sqrtau * t ** (tmp + 0.5_dp) * exp(-t) * x

    end if cgf
  end function rcsv_gamma

  recursive function rcsv_ln_gamma(z) result(ln_gamma)
    ! Calculates the value of the natural logarithm of the gamma function for gamma = 8 and N = 9
    implicit none
    complex(dpc), intent(in) :: z
    complex(dpc) :: ln_gamma
    complex(dpc) :: tmp, t
    real(dp) :: x
    integer  :: i
    real(dp), parameter, dimension(9) :: coef = [0.99999999999980993_dp, 676.5203681218851_dp, -1259.1392167224028_dp, &
    771.32342877765313_dp, -176.61502916214059_dp, 12.507343278686905_dp,&
    -0.13857109526572012_dp, 9.9843695780195716d-6, 1.5056327351493116d-7]
    real(dp), parameter :: lnsqrtau = log(sqrt(tau)), lnpi = log(pi)

    tmp = z

    ! Calculate Ln(Gamma), it takes the main branch cut.
    clg: if (real(tmp) < 0.5_dp) then

      ln_gamma = lnpi - log(sin(pi*tmp)) - rcsv_ln_gamma(1._dp-tmp)

    else clg

       tmp = tmp - 1._dp
       x   = coef(1)

       ! Calculate tge Arithmetic Progression and Divide.
       apd: do i = 2, size(coef)
         x = x + coef(i) / (tmp + i - 1)
       end do apd

       t = tmp + size(coef) - 1.5_dp

       ln_gamma = lnsqrtau + log(t)*(tmp + 0.5_dp) -t + log(x)

    end if clg

  end function rcsv_ln_gamma
end module nrtools
