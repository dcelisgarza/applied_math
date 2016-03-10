program solar_system_main
  use orbital_mech
  implicit none
  real(dp) :: kep_el(20,10) ! Keplerian elements of the bodies of the solar system.
  type(orb_par) :: orbits(10) ! Orbits of the solar system.
  character(:), allocatable :: filename
  integer  :: i
  real(dp) :: t, dt, finish ! time, increment in time.

  dt = 1._dp
  finish =  690._dp


  ! Read keplerian elements.
  !filename = 'kep_el.txt'
  !call read_kep_el(filename, size(orbits), 1, kep_el)

  ! Evolve keplerian elements per century.
  !call kep_el_teph(kep_el, orbits, 2451545._dp)

  ! Calculate the state vectors.
  !call state_vec(orbits,1d-6)

  ! Evolve the trajectories with time.
  call integrate_orbits(orbits, 0._dp, finish, dt)
end program solar_system_main
