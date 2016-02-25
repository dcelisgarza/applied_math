module orbital_mech
  use nrtype
  use find_roots

  type orb_par
    real(dp) :: a  ! Major semi-axis.
    real(dp) :: ea ! Eccentric anomaly.
    real(dp) :: ec ! Eccentricity.
    real(dp) :: ma ! Mean anomaly.
    real(dp) :: mu ! Gravitational parameter.
    real(dp) :: ta ! True anomaly.
    real(dp) :: dt ! Time step.
    real(dp) :: n  ! Mean angular motion := sqrt(mu/a**3)
    real(dp) :: r  ! Distance from the focus.
    real(dp) :: o(4) ! Orbital position and velocity.
    real(dp) :: x(6) ! Cartesian position and velocity.
    real(dp) :: eu(3)! Euler angles. eu(1) = Omega, eu(2) = i, eu(3) = omega
    real(dp) :: reu(6) ! Euler angles in the rotation matrix.
  end type orb_par

contains
  subroutine eorb(orbits,nr_tol)
    implicit none
    ! Evolves the orbital path.
    type(orb_par), intent(inout) :: orbits(:)
    real(dp), intent(in)         :: nr_tol
    real(dp), dimension(size(orbits)) :: c1, c2, cosea, eao2
    integer :: i

    ! Evolve Mean Anomaly along the time step.
    orbits % ma = mod(orbits % ma + orbits % dt * orbits % n, tau)
    ! Solve Kepler's equation for all orbits to calculate the eccentric anomaly.
    call orb_new_rap(orbits,nr_tol)

    ! Calculate true anomaly.
    eao2  = orbits % ea * .5_dp
    cosea = cos( orbits % ea )
    orbits % ta = 2._dp * atan2( sqrt(1-orbits%ec)*cos( eao2 ), &
    sqrt(1 + orbits % ec) * sin( eao2 ) )

    ! Calculate distance from the focus.
    orbits % r = orbits % a * (1._dp - orbits % ec * cosea )
    ! Calculating orbital positions.
    ! Calculating x in orbital coordinates.
    orbits % o(1) = orbits % r * cos( orbits % ta )
    ! Calculating y in orbital coordiantes.
    orbits % o(2) = orbits % r * sin( orbits % ta )

    ! Calculating orbital velocities.
    c1 = sqrt( orbits % mu * orbits % a ) / orbits % r
    c2 = sqrt( 1 - orbits % ec * orbits % ec )
    ! Calculating dx/dea.
    orbits % o(3) = - c1 * sin( orbits % ea )
    ! Calculating dy/dea.
    orbits % o(4) = c1 * c2 * cosea

    orb_to_cart: do i = 1, size(orbits)
      orbits(i) % x(1:3) = orbits(i) % o(1) * orbits(i) % reu(1:3) + orbits(i) % o(2) * orbits(i) % reu(4:6)
      orbits(i) % x(4:6) = orbits(i) % o(3) * orbits(i) % reu(1:3) + orbits(i) % o(4) * orbits(i) % reu(4:6)
    end do orb_to_cart
  end subroutine eorb

  subroutine kepler_eqn(kepeq,dkepeq,orbits)
    implicit none
    type(orb_par), intent(in)                      :: orbits(:)
    real(dp), intent(out), dimension(size(orbits)) :: kepeq, dkepeq

    ! We need to know when this is zero. From here we obtain the eccentric anomaly.
    kepeq  = orbits % ea - orbits % ec * sin(orbits % ea) - orbits % ma
    dkepeq = 1 - orbits % ec * cos(orbits % ea)
  end subroutine kepler_eqn

  subroutine orb_new_rap(orbits,tol)
    implicit none
    type(orb_par), intent(inout)      :: orbits(:)
    real(dp), dimension(size(orbits)) :: kepeq, dkepeq
    real(dp), intent(in)              :: tol
    real(dp) :: ea_1(size(orbits))
    integer  :: obj

    obj  = size(orbits)

    new_rap: do
      ea_1 = orbits % ea
      call kepler_eqn(kepeq,dkepeq,orbits)
      orbits % ea = orbits % ea - kepeq/dkepeq
      if ( 0.5_dp * abs( sum(orbits % ea)/obj - sum(ea_1)/obj ) < tol) return
    end do new_rap
  end subroutine orb_new_rap
end module orbital_mech
