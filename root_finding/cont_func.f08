module continuous_functions
  use numbers
contains
  pure function polynomial(x)
    implicit none
    real(dp), intent(in),optional :: x
    real(dp)             :: polynomial

    polynomial = x**3. - x - 2.
    !(x-1.)*(x+1.)*(x-2.)*(x+2.)*(x-3.)*(x+3.)
    !1.*x**6.-1./3.*x**5.+ 1./5.*x**4.-1./7.*x**3.+1./9.*x**2.-1./11.*x+1./13.
  end function polynomial
end module continuous_functions
