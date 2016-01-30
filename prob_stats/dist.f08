module distributions
  use numbers
contains
  pure elemental function radial_sine_dist(r)
    implicit none
    real(dp), intent(in) :: r
    real(dp) radial_sine_dist
    radial_sine_dist = sin(r)**2./r**3.
  end function radial_sine_dist

end module distributions
