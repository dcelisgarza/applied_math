module orbital_mech
  use nrtype
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
    real(dp) :: eu(3)! Euler angles. eu(1) = Omega, eu(2) = i, eu(3) = omega
    real(dp) :: reu(6) ! Euler angles in the rotation matrix.
  end type orb_par

contains
  subroutine read_kep_el(filename, bodies, kep_el)
    implicit none
    ! Reads Keplerian Elements.
    character(:), allocatable, intent(in) :: filename
    integer, intent(in)      :: bodies
    real(dp), intent(out)    :: kep_el(:,:)
    integer                  :: open_status, close_status, i, j

    open(unit = 1, file = filename, status = 'old', iostat = open_status, action = 'read', position = 'rewind')
    if ( open_status /= 0 ) then
      print *, 'Could not open ',filename,' for reading.', &
      'unit = ', 1
      stop
    endif

    read_file1: do i = 1, bodies
      read(1,*) kep_el(1:6,i)
      read(1,*) kep_el(7:12,i)
    end do read_file1

  end subroutine read_kep_el

  subroutine kep_el_teph(kep_el, orbits, t_eph)
    implicit none
    ! Calculates the Keplerian Elements at a given Julian Epoch.
    real(dp), intent(inout)    :: kep_el(:,:)
    real(dp), intent(in)       :: t_eph ! Julian Ephemeris Date
    type(orb_par), intent(out) :: orbits(:)
    real(dp)                   :: t

    ! Calculating T for evolving orbital elements in time.
    t = (t_eph-2451545._dp)/36525

    ! Evolve orbital elements in time.
    kep_el(1:6,:)   = kep_el(1:6,:)   + kep_el(7:12,:)  * t
    kep_el(13:16,:) = kep_el(13:16,:) + kep_el(17:20,:) * t

    ! Calculate orbital parameters.
    orbits % a     = kep_el(1,:)
    orbits % ec    = kep_el(2,:)
    orbits % eu(1) = kep_el(6,:) * pi/180._dp
    orbits % eu(2) = kep_el(3,:) * pi/180._dp
    orbits % eu(3) = (kep_el(5,:) - kep_el(6,:)) * pi/180._dp
    orbits % ma    = kep_el(4,:) - kep_el(5,:)        &
    + kep_el(13,:)*t*t                 &
    + kep_el(14,:)*cos(kep_el(16,:)*t) &
    + kep_el(15,:)*sin(kep_el(16,:)*t)
    orbits % ma    = mod(orbits % ma, tau) - pi

  end subroutine kep_el_teph

  subroutine state_vec(orbits,nr_tol)
    implicit none
    ! Calculates position and velocity state vectors.
    type(orb_par), intent(inout)      :: orbits(:)
    real(dp), intent(in)              :: nr_tol
    real(dp), dimension(size(orbits)) :: c1, c2, cosea, eao2
    integer :: i

    ! Solve Kepler's equation for all orbits to calculate the eccentric anomaly.
    call orb_new_rap(orbits,nr_tol)

    ! Calculate true anomaly.
    eao2  = orbits % ea * .5_dp
    cosea = cos( orbits % ea )
    orbits % ta = 2._dp * atan2( sqrt( 1 - orbits % ec )*cos( eao2 ), &
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

    ! Call the rotation matrix for the euler angles.
    call rotate_orb_cart(orbits)

    orb_to_cart: do i = 1, size(orbits)
      orbits(i) % x(1:3) = orbits(i) % o(1) * orbits(i) % reu(1:3) + orbits(i) % o(2) * orbits(i) % reu(4:6)
      orbits(i) % x(4:6) = orbits(i) % o(3) * orbits(i) % reu(1:3) + orbits(i) % o(4) * orbits(i) % reu(4:6)
    end do orb_to_cart
  end subroutine state_vec

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

  subroutine rotate_orb_cart(orbits)
    implicit none
    type(orb_par), intent(inout) :: orbits(:)
    orbits % reu(1) =  cos(orbits % eu(1)) * cos(orbits % eu(3)) - sin(orbits % eu(1)) * cos(orbits % eu(2)) * sin(orbits % eu(3))
    orbits % reu(2) =  sin(orbits % eu(1)) * cos(orbits % eu(3)) + cos(orbits % eu(1)) * cos(orbits % eu(2)) * sin(orbits % eu(3))
    orbits % reu(3) =  sin(orbits % eu(2)) * sin(orbits % eu(3))
    orbits % reu(4) = -cos(orbits % eu(1)) * sin(orbits % eu(3)) - sin(orbits % eu(1)) * cos(orbits % eu(2)) * cos(orbits % eu(3))
    orbits % reu(6) =  sin(orbits % eu(2)) * sin(orbits % eu(1))
    orbits % reu(5) = -sin(orbits % eu(1)) * sin(orbits % eu(3)) + cos(orbits % eu(1)) * cos(orbits % eu(2)) * cos(orbits % eu(3))
  end subroutine rotate_orb_cart
end module orbital_mech
