module distributions
  use nrtype
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

  pure function radial_exp_decay_cos(radius,angle)
    implicit none
    real(dp), intent(in), optional :: radius
    real(dp), intent(in), optional :: angle
    real(dp) radial_exp_decay_cos
    radial_exp_decay_cos = cos(radius)**2.*exp(-radius)
  end function radial_exp_decay_cos

  pure function bivariate_gaussian(x, aux)
    implicit none
    real(dp), intent(in), optional :: x(:), aux(:)

    real(dp) :: bivariate_gaussian
    real(dp) :: sig_x_y, rho_fact, d_x_mu_x, d_y_mu_y

    ! Equivalence.
    ! x = x(1)
    ! y = x(2)
    ! rho = aux(1)
    ! mu_x = aux(2)
    ! mu_y = aux(3)
    ! sig_x = aux(4)
    ! sig_y = aux(5)

    sig_x_y = 1./aux(4)*aux(5)
    rho_fact = 1./(1-aux(1)**2.)
    d_x_mu_x = x(1)-aux(2)
    d_y_mu_y = x(2)-aux(3)

    bivariate_gaussian = 1./tau * sig_x_y * sqrt(rho_fact) * &
    exp( &
    - 0.5 * rho_fact * &
    ( &
      d_x_mu_x ** 2. / aux(4) ** 2. &
    + d_y_mu_y ** 2. / aux(5) ** 2. &
    - 2. * aux(1) * d_x_mu_x * d_y_mu_y / ( aux(4) * aux(5) ) &
    ) &
    )
  end function bivariate_gaussian
end module distributions
