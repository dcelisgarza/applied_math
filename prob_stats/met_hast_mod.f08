module metropolis_hastings_mod
  use nrtype
  use distributions
contains

  subroutine n_dim_met_hast(dist,n_dim,aux,mutagen,auto_corr,nbr_samples,sample,file_unit)

    implicit none
    interface n_dim_distribution
      function dist(x,aux)
        use nrtype
        implicit none
        real(dp), intent(in), optional :: x(:)
        real(dp), intent(in), optional :: aux(:)
        real(dp)                       :: dist
      end function dist
    end interface n_dim_distribution
    integer, intent(in) :: n_dim
    real(dp), intent(in) :: aux(:)
    real(dp), intent(inout), optional :: sample(n_dim) ! Sample.
    integer, intent(in) :: auto_corr, nbr_samples ! Auto correlation, Number of samples, file unit.
    integer, intent(in), optional :: file_unit
    integer  :: i, j      ! counters
    real(dp) :: mutagen(:)
    real(dp) :: accept ! Selection criterion
    real(dp) :: current_sample(n_dim) ! Current sample.
    real(dp) :: x(2,n_dim), dx(n_dim)
    real(dp) :: ZBQLU01, ZBQLUAB
    external ZBQLU01, ZBQLUAB

    call ZBQLINI(0)
    ! Not a self-starter.
    current_sample = 0.
    i = 0

    collect_samples: do while (i < auto_corr * nbr_samples)
      ! Add a flag in the parameters so that you can change the mutation function depending on what type of distribution you want to sample. You may also need to change or add the loops.
      ! Dimension differentials.
      dim_diff: do j = 1, n_dim
        dx(j) = ZBQLUAB(mutagen(1),mutagen(2))
      end do dim_diff

      ! Current point.
      x(1,:) = current_sample

      ! Next point.
      x(2,:) = x(1,:) + dx

      ! Check whether our new coordinates fit the distribution.
      ! Min( f(X_2)/f(X_1), 1 )
      accept = min(dist(x(2,:),aux) / dist(x(1,:),aux), 1.)
      acceptance: if (accept > ZBQLU01(0.)) then
        current_sample = x(2,:)
        i = i + 1
        take_sample: if (mod(i, auto_corr) == 0) then
          write_file: if (present(file_unit)) then
            write(file_unit, *) current_sample
          else write_file
            sample = current_sample
          end if write_file
        end if take_sample
      end if acceptance
    end do collect_samples

    if (present(file_unit)) close(file_unit, status = 'keep')

  end subroutine n_dim_met_hast

  subroutine circular_met_hast(dist,radial,angular,auto_corr,nbr_samples,sample,file_unit)

    implicit none
    interface polar_distribution
      function dist(radius,angle)
        use nrtype
        implicit none
        real(dp), intent(in), optional :: radius
        real(dp), intent(in), optional :: angle
        real(dp)                       :: dist
      end function dist
    end interface polar_distribution
    real(dp), intent(in) :: radial, angular(2) ! Radial random number, Angular random number.
    real(dp), intent(inout), optional :: sample(2) ! Sample.
    integer, intent(in) :: auto_corr, nbr_samples ! Auto correlation, Number of samples, file unit.
    integer, intent(in), optional :: file_unit
    integer  :: i      ! counter
    real(dp) :: accept ! Selection criterion
    real(dp) :: current_sample(2) ! Current sample.
    real(dp) :: r(2), th(2), x(2), y(2) ! r, x, y.
    real(dp) :: radius, theta
    real(dp) :: ZBQLU01, ZBQLEXP, ZBQLUAB
    external ZBQLU01, ZBQLEXP, ZBQLUAB

    call ZBQLINI(0)
    ! Not a self-starter.
    current_sample = 0.
    i = 0

    collect_samples: do while (i < auto_corr * nbr_samples)
      ! Random radius and angle.
      radius = ZBQLEXP(radial)
      theta  = ZBQLUAB(angular(1),angular(2))

      ! Current point.
      x(1) = current_sample(1)
      y(1) = current_sample(2)

      ! Next point.
      x(2) = x(1) + radius * cos(theta)
      y(2) = y(1) + radius * sin(theta)

      ! Calculating polar coordinates for distribution.
      r  = sqrt(x**2. + y**2.)
      th = atan2(y, x)

      ! Check whether our new coordinates fit the distribution.
      ! Min( f(X_2)/f(X_1), 1 )
      accept = min( dist(r(2),th(2)) / dist(r(1),th(1)), 1. )
      acceptance: if (accept > ZBQLU01(0.)) then
        current_sample(1) = x(2)
        current_sample(2) = y(2)
        i = i + 1
        take_sample: if (mod(i, auto_corr) == 0) then
          write_file: if (present(file_unit)) then
            write(file_unit, *) current_sample
          else write_file
            sample = current_sample
          end if write_file
        end if take_sample
      end if acceptance
    end do collect_samples

    if (present(file_unit)) close(file_unit, status = 'keep')
  end subroutine circular_met_hast

  !=======================================================================================!

  subroutine polar_met_hast(dist,radial,angular,auto_corr,nbr_samples,sample,file_unit)

    implicit none
    interface polar_distribution
      function dist(radius,angle)
        use nrtype
        implicit none
        real(dp), intent(in), optional :: radius
        real(dp), intent(in), optional :: angle
        real(dp)                       :: dist
      end function dist
    end interface polar_distribution
    real(dp), intent(in) :: radial, angular(2) ! Radial random number, Angular random number.
    real(dp), intent(inout), optional :: sample(2) ! Sample.
    integer, intent(in) :: auto_corr, nbr_samples ! Auto correlation, Number of samples, file unit.
    integer, intent(in), optional :: file_unit
    integer  :: i      ! counter
    real(dp) :: accept ! Selection criterion
    real(dp) :: current_sample(2) ! Current sample.
    real(dp) :: r(2), th(2) ! radius and theta
    real(dp) :: radius, theta
    real(dp) :: ZBQLU01, ZBQLEXP, ZBQLUAB
    external ZBQLU01, ZBQLEXP, ZBQLUAB

    call ZBQLINI(0)
    ! Not a self-starter.
    current_sample = 0.
    i = 0

    collect_samples: do while (i < auto_corr * nbr_samples)
      ! Random radius and angle.
      radius = ZBQLEXP(radial)
      theta  = ZBQLUAB(angular(1),angular(2))

      ! Current point.
      r(1)  = current_sample(1)
      th(1) = current_sample(2)

      ! Next point.
      r(2)  = r(1)  + radius
      th(2) = th(1) + theta

      ! Check whether our new coordinates fit the distribution.
      ! Min( f(X_2)/f(X_1), 1 )
      accept = min( dist(r(2),th(2)) / dist(r(1),th(1)), 1. )
      acceptance: if (accept > ZBQLU01(0.)) then
        current_sample(1) = r(2)
        current_sample(2) = th(2)
        i = i + 1
        take_sample: if (mod(i, auto_corr) == 0) then
          write_file: if (present(file_unit)) then
            write(file_unit, *) current_sample
          else write_file
            sample = current_sample
          end if write_file
        end if take_sample
      end if acceptance
    end do collect_samples

    if (present(file_unit)) close(file_unit, status = 'keep')

  end subroutine polar_met_hast

  !=======================================================================================!

  subroutine spherical_met_hast(dist,radial,angular,auto_corr,nbr_samples,sample,file_unit)

    implicit none
    interface spherical_distribution
      function dist(radius,theta,phi)
        use nrtype
        implicit none
        real(dp), intent(in), optional :: radius
        real(dp), intent(in), optional :: theta
        real(dp), intent(in), optional :: phi
        real(dp)                       :: dist
      end function dist
    end interface spherical_distribution
    real(dp), intent(in) :: radial, angular(4) ! Radial random number, Angular random numbers.
    real(dp), intent(inout), optional :: sample(3) ! Sample.
    integer, intent(in) :: auto_corr, nbr_samples ! Auto correlation, Number of samples, file unit.
    integer, intent(in), optional :: file_unit
    integer  :: i      ! counter
    real(dp) :: accept ! Selection criterion
    real(dp) :: current_sample(3) ! Current sample.
    real(dp) :: x(2), y(2), z(2), r(2), th(2), ph(2) ! x, y, z, radius, theta, phi.
    real(dp) :: radius, theta, phi
    real(dp) :: ZBQLU01, ZBQLEXP, ZBQLUAB
    external ZBQLU01, ZBQLEXP, ZBQLUAB

    call ZBQLINI(0)
    ! Not a self-starter.
    current_sample = 0.
    i = 0

    collect_samples: do while (i < auto_corr * nbr_samples)
      ! Random radius and angle.
      radius = ZBQLEXP(radial)
      theta  = ZBQLUAB(angular(1),angular(2)) ! theta \in [0, 2pi]
      phi    = ZBQLUAB(angular(3),angular(4)) ! phi \in [0, pi]

      ! Current point.
      x(1) = current_sample(1)
      y(1) = current_sample(2)
      z(1) = current_sample(3)

      ! Next point.
      x(2) = x(1) + radius * cos(theta) * sin(phi)
      y(2) = y(1) + radius * sin(theta) * sin(phi)
      z(2) = z(1) + radius * cos(phi)

      ! Spherical coordinates for the distributions.
      r  = sqrt(x**2. + y**2. + z**2.)
      th = atan2(y, x)
      ph = acos(z/r)

      ! Check whether our new coordinates fit the distribution.
      ! Min( f(X_2)/f(X_1), 1 )
      accept = min( dist(r(2),th(2),ph(2)) / dist(r(1),th(1),ph(1)), 1. )
      acceptance: if (accept > ZBQLU01(0.)) then
        current_sample(1) = r(2)
        current_sample(2) = th(2)
        i = i + 1
        take_sample: if (mod(i, auto_corr) == 0) then
          write_file: if (present(file_unit)) then
            write(file_unit, *) current_sample
          else write_file
            sample = current_sample
          end if write_file
        end if take_sample
      end if acceptance
    end do collect_samples

    if (present(file_unit)) close(file_unit, status = 'keep')

  end subroutine spherical_met_hast

  !=======================================================================================!

end module metropolis_hastings_mod
