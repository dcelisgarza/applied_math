module orbital_mech
  use nrtype
  use ode_int, only : velocity_verlet

  type orb_par
    real(dp) :: a  ! Major semi-axis.
    real(dp) :: ea ! Eccentric anomaly.
    real(dp) :: ec ! Eccentricity.
    real(dp) :: ma ! Mean anomaly.
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
  subroutine integrate_orbits(orbits, start, end, dt)
    implicit none
    real(dp), intent(in) :: start, end, dt
    type(orb_par), intent(inout) :: orbits(:)
    real(dp),dimension(6*size(orbits)) :: xi, xf ! initial orbital positions and velocities
    real(dp) :: t  ! time
    integer  :: i, nbod, n

    nbod = size(orbits)

    open(1, file = 'vectors.txt')
    read_init: do i = 1, nbod
      read(1,*) orbits(i) % x(1:3) ! Positions, X, Y, Z.
      read(1,*) orbits(i) % x(4:6) ! Velocities Vx, Vy, Vz.
    end do read_init
    close(1)

    open(unit = 1, file = 'solar_system.dat', status = 'replace', action = 'write', position = 'rewind')

    n = 0
    ! Map orbital positions to a 1D array with 6n entries, where n = number of bodies.
    do i = 1, nbod
      ! xi(1:3*nbod) = x, y, z.
      xi(i+n:i+n+2) = orbits(i)%x(1:3)
      ! xi(3*nbod+1:6*nbod) = vx, vy, vz.
      xi(i+n+nbod*3:i+n+nbod*3+2) = orbits(i)%x(4:6)
      n = n + 2
    end do

      write(1,*) xi

    do while( t < end)
        call velocity_verlet(accel_solar_system,t,xi,xf,dt)
        t = t + dt
        if (mod(t, 2.*dt) == 0.) write(1,*) xf
        xi = xf
    end do
  end subroutine integrate_orbits

  subroutine read_kep_el(filename, bodies, date_flag, kep_el)
    implicit none
    ! Reads Keplerian Elements.
    character(:), allocatable, intent(in) :: filename
    integer, intent(in)      :: bodies, date_flag
    real(dp), intent(out)    :: kep_el(:,:)
    integer                  :: open_status, close_status, i, j

    open(unit = 1, file = filename, status = 'old', iostat = open_status, action = 'read', position = 'rewind')
    if ( open_status /= 0 ) then
      print *, 'Could not open ',filename,' for reading.', &
      'unit = ', 1
      stop
    endif

    epoch: if (date_flag == 1) then
      read_file1: do i = 1, bodies
        read(1,*) kep_el(1:10,i)
        read(1,*) kep_el(11:20,i)
      end do read_file1
    else epoch
      do j = 1, 19
        read(1,*)
      end do
      read_file2: do i = 1, bodies
        read(1,*) kep_el(1:10,i)
        read(1,*) kep_el(11:20,i)
      end do read_file2
    end if epoch
    close(unit = 1)
  end subroutine read_kep_el

  subroutine kep_el_teph(kep_el, orbits, t_eph)
    implicit none
    ! Calculates the Keplerian Elements at a given Julian Epoch.
    real(dp), intent(inout)    :: kep_el(:,:)
    real(dp), intent(in)       :: t_eph ! Julian Ephemeris Date
    type(orb_par), intent(out) :: orbits(:)
    real(dp)                   :: t

    ! Calculating T for evolving orbital elements in time.
    t = (t_eph-2451545._dp)/36525._dp

    ! Evolve orbital elements in time.
    kep_el(1:6,:)   = kep_el(1:6,:)   + kep_el(7:12,:)  * t
    kep_el(13:16,:) = kep_el(13:16,:) + kep_el(17:20,:) * t

    ! Calculate orbital parameters.
    orbits % a     = kep_el(1,:)
    orbits % ec    = kep_el(2,:)
    orbits % eu(1) = kep_el(6,:) * pi/180._dp
    orbits % eu(2) = kep_el(3,:) * pi/180._dp
    orbits % eu(3) = (kep_el(5,:) - kep_el(6,:)) * pi/180._dp
    orbits % ma    = kep_el(4,:)  - kep_el(5,:)        &
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
    real(dp), parameter :: mu = 1.32712440018d9/1.49597870700_dp, secoday = 1._dp/86400._dp ! seconds
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
    c1 = sqrt( mu * orbits % a ) / orbits % r
    c2 = sqrt( 1 - orbits % ec * orbits % ec )
    ! Calculating dx/dea.
    orbits % o(3) = - c1 * sin( orbits % ea )*secoday
    ! Calculating dy/dea.
    orbits % o(4) = c1 * c2 * cosea*secoday

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
    ! Initial value of the eccentric anomaly is the mean anomaly.
    orbits % ea = orbits % ma
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

  subroutine accel_solar_system(x,y,dydx)
    implicit none
    real(dp), intent(in)  :: x, y(:)
    real(dp), intent(out) :: dydx(:)
    real(dp), parameter   :: kmau   = 1.495978707E+8, &   ! km/Au
                             km3au3  = kmau*kmau*kmau, &  ! km^3/Au^3
                             sd  = 8.64E+4, &             ! s/d
                             s2d2 = sd*sd                 ! s^2/d^2
    real(dp), parameter, dimension(10) :: g = [1.3271244004193938E+11, 2.203209E+4, 3.2485863E+5, &
                                               3.9860044E+5, 4.28283E+4, 1.26686511E+8, 3.79312078E+7,&
                                               5.793966E+6, 6.835107E+6, 872.4]/km3au3*s2d2
    integer  :: i, j, m, n ! i and j are coutners. m and n are counters that increase by two
    real(dp) :: norm, vec(3)
    integer  :: nbod

    nbod = size(y)/6
    ! Initial acceleration is zero because we need to add the acceleration caused by all interactions.
    dydx = 0._dp

    ! M increases by two every iteration of the loop for the planets. This is due to mapping an n x 6 matrix to a vector with 6n entries, where n = # of bodies.
    n = 0
    loop_cel_bodies: do i = 1, nbod
      ! N inreases by two every iteration of the acceleration. This is due to mapping an n x 6 matrix to a vector with 6n entries, where n = # of bodies.
      m = 0
      loop_accel: do j = 1, nbod
        ! Check for self-interactions, and skip them. We need to add two to m because the potential energy difference between a celestial body and itself is zero, therefore the acceleration with respect to itself is zero.
        if (j == i) then
          m = m + 2
          cycle loop_accel
        end if
        ! Calculate the components of the direction vector.
        vec  =  y(j + m : j + 2 + m) - y(i + n: i + 2 + n)
        ! Calculate the norm of the direction vector.
        norm = sqrt( sum(vec*vec) )
        ! Cube the norm.
        norm = norm*norm*norm
        ! Calculate the the total acceleration for a celestial body.
        dydx(i + n: i + 2 + n) = g(j)*vec / norm + dydx(i + n: i + 2 + n)
        m = m + 2
      end do loop_accel
      n = n + 2
    end do loop_cel_bodies
  end subroutine accel_solar_system
end module orbital_mech
