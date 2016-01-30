module distributions
  use numbers
contains
  pure function radial_sine_dist(radius,angle)
    implicit none
    real(dp), intent(in), optional :: radius
    real(dp), intent(in), optional :: angle
    real(dp) radial_sine_dist
    radial_sine_dist = sin(radius)**2./radius**3.
  end function radial_sine_dist

  pure function radial_sine_cos_dist(radius,angle)
    implicit none
    real(dp), intent(in), optional :: radius
    real(dp), intent(in), optional :: angle
    real(dp) radial_sine_cos_dist
    radial_sine_cos_dist = sin(radius)**2.*cos(radius)**2./radius**3.
  end function radial_sine_cos_dist

  pure function exp_decay_cos(radius,angle)
    implicit none
    real(dp), intent(in), optional :: radius
    real(dp), intent(in), optional :: angle
    real(dp) exp_decay_cos
    exp_decay_cos = cos(radius)**2.*exp(-radius)
  end function exp_decay_cos

end module distributions
