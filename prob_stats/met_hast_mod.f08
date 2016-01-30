module metropolis_hastings_mod
  use numbers
  use distributions
contains

  subroutine circular_met_hast(dist,radial,angular,auto_corr,nbr_samples,sample,file_unit)

    implicit none
    interface radial_distribution
      function dist(radius,angle)
        use numbers
        implicit none
        real(dp), intent(in), optional :: radius
        real(dp), intent(in), optional :: angle
        real(dp)                       :: dist
      end function dist
    end interface radial_distribution
    real(dp), intent(in) :: radial, angular(2) ! Radial random number, Angular random number.
    real(dp), intent(inout), optional :: sample(2) ! Sample.
    integer, intent(in) :: auto_corr, nbr_samples ! Auto correlation, Number of samples, file unit.
    integer, intent(in), optional :: file_unit
    integer :: i      ! counter
    real(dp) :: accept ! Selection criterion
    real(dp) :: current_sample(2) ! Current sample.
    real(dp) :: r(2), x(2), y(2) ! x, y.
    real(dp) :: distance, direction, radius(2), angle(2)
    real(dp) :: ZBQLU01, ZBQLEXP, ZBQLUAB
    external ZBQLU01, ZBQLEXP, ZBQLUAB

    call ZBQLINI(0)
    ! Not a self-starter.
    current_sample = 0.
    i = 0

    collect_samples: do while (i < auto_corr * nbr_samples)
      ! Random radius and angle.
      distance  = ZBQLEXP(radial)
      direction = ZBQLUAB(angular(1),angular(2))

      ! Current point.
      x(1) = current_sample(1)
      y(1) = current_sample(2)

      ! Next point.
      x(2) = x(1) + distance*cos(direction)
      y(2) = y(1) + distance*sin(direction)

      ! Check whether our new coordinates fit the distribution.
      ! Min( f(X_2)/f(X_1), 1 )
      radius = sqrt(x**2. + y**2.)
      angle  = atan2(y, x)
      accept = min( dist(radius(2),angle(2)) / dist(radius(1),angle(1)), 1. )
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

end module metropolis_hastings_mod
