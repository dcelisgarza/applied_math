program main
  use nrtype
  use distributions

  implicit none
  integer  :: i, ac, nbr_samples   ! i = counter, ac = autocorrelation, number of samples.
  real(dp) :: accept               ! accept = acceptance criterion.
  real(dp) :: direction, distance
  real(dp) :: current_sample(2)    ! current sample, could or could not be rejected.
  real(dp) :: r(2), x(2), y(2)     ! values of r, x and y.

  real(dp) :: ZBQLUAB, ZBQLEXP, ZBQLU01
  external ZBQLINI, ZBQLUAB, ZBQLEXP, ZBQLU01

  call ZBQLINI(0)

  current_sample = 0.!1d-8

  i  = 0
  ac = 5
  nbr_samples = 1e6
  open(unit = 1, file = 'data.dat')

  collect_samples: do while (i < ac * nbr_samples)
    ! Direction and distance.
    direction = ZBQLUAB(0.,tau)
    distance  = -DLOG(ZBQLU01(0.0D0))*2.!ZBQLEXP(2.)

    ! Current point.
    x(1) = current_sample(1)
    y(1) = current_sample(2)

    ! Next point.
    x(2) = x(1) + distance*cos(direction)
    y(2) = y(1) + distance*sin(direction)

    ! Check whether it fits the distribution.
    ! Min( f(X_2)/f(X_1), 1 )
    r = sqrt(x**2. + y**2.)
    accept = min( radial_sine_dist(r(2)) / radial_sine_dist(r(1)), 1.)
    acceptance: if (accept > ZBQLU01(0.)) then
      current_sample(1) = x(2)
      current_sample(2) = y(2)
      i = i + 1
      take_sample: if (mod(i, ac) == 0) then
        write(1,*) current_sample
      end if take_sample
    end if acceptance
  end do collect_samples
  close (1, status='keep')

  call system( "gnuplot plot_met_hast.plt" )
  !call system( "" )

end program main
