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
    integer(i16)        :: bicoef

    ! Identify Key Scenarios.
    iks: if (k > n) then
      write(*,*) " Error: nrtools: rcsv_bicoef: iks: n must be greater than k. k = ", k, ", n = ", n
      return
    else if (k == 0) then iks
      bicoef = 1
    else if (k > n/2) then iks
      ! Recurse function so as to reduce the number of operations to at most n/2 at the numerator and denominator.
      bicoef = rcsv_bicoef(n, n-k)
    else iks
      ! Calculate binomial coefficient.
      bicoef = n/k * rcsv_bicoef(n-1, k-1)
    end if iks

  end function rcsv_bicoef
end module nrtools
