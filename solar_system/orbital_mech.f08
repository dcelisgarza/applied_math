module orbital_mech
  use nrtype
  use find_roots
  use ode_int

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
    real(dp) :: eu(3)! Euler angles.
  end type orb_par

contains
  subroutine eorb(orbits,tol)
    implicit none
    ! Evolves the orbital path.
    type(orb_par), intent(inout) :: orbits(:)
    real(dp), intent(in)         :: tol

    ! Evolve Mean Anomaly along the time step.
    orbits % ma = orbits % ma + orbits % dt * orbits % n

    ! Solve Kepler's equation for all orbits.
     call newt_rap(orbits,tol)

     ! Calculate true anomaly.
     orbits % ta = 2._dp * atan2( sqrt(1-orbits%ec)*cos(orbits%ea*.5_dp), &
                                  sqrt(1+orbits%ec)*sin(orbits%ea*.5_dp) )

     ! Calculate distance from the focus.
     orbits % r = orbits % a * (1._dp - orbits % ec * cos(orbits % ea) )
     ! Calculating orbital positions.
     ! Calculating x in orbital coordinates.
     orbits % o(1) = orbits % r * cos(orbits % ta)
     ! Calculating y in orbital coordiantes.
     orbits % o(2) = orbits % r * sin(orbits % ta)

     ! Calculating orbital velocities.
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
    type(orb_par), intent(inout)      :: orbits
    real(dp), dimension(size(orbits)) :: kepeq, dkepeq
    real(dp), intent(in)              :: err
    interface kepler
      subroutine kepler_eqn(kepeq,dkepeq,orbits)
        implicit none
        type(orb_par), intent(in)                      :: orbits(:)
        real(dp), intent(out), dimension(size(orbits)) :: kepeq, dkepeq
      end subroutine kepler_eqn
    end interface kepler
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
