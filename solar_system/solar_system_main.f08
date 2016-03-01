program solar_system_main
  use orbital_mech
  implicit none
  real(dp) :: kep_el(20,10) ! Keplerian elements of the bodies of the solar system.
  type(orb_par) :: orbits(10) ! Orbits of the solar system.
  character(:), allocatable :: filename
  integer  :: i
  real(dp) :: t, dt ! time, increment in time.

  dt = 1._dp


  ! Read keplerian elements.
  filename = 'kep_el.txt'
  call read_kep_el(filename, size(orbits), 1, kep_el)
  ! Evolve keplerian elements per century.
  call kep_el_teph(kep_el, orbits, 2451545._dp)
  ! Calculate the state vectors.
  call state_vec(orbits,1d-6)
  ! Evolve the trajectories with time.
!  open(unit = 1, file = 'solar_system.dat', status = 'replace', action = 'write', position = 'rewind')
!  do i = 1, 10
!    write(1,100,advance='no') orbits(i) % x(1:3)
!  end do

  orbits(1)%x=0._dp
  !print*, orbits(2)%x, '|', orbits(3)%x
  !print*, ''
  !print*, orbits(3)%x - orbits(2)%x
  !print*, '------------'
  call integrate_orbits(orbits, 0._dp, 1._dp, 0.1_dp)

  100 format (10(3F8.4))

  !  do i = 1, 10
  !    write(*,*) kep_el(1:10,i)
  !    write(*,*) kep_el(11:20,i)
  !  end do
end program solar_system_main
