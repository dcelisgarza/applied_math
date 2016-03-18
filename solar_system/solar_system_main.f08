program solar_system_main
  use orbital_mech
  use plot
  implicit none
  real(dp) :: kep_el(20,10) ! Keplerian elements of the bodies of the solar system.
  type(orb_par) :: orbits(10) ! Orbits of the solar system.
  character(:), allocatable :: filename
  integer  :: i
  real(dp) :: t, dt, start, finish, rec ! time, increment in time.

  dt = 1._dp
  start = 0._dp
  finish =  690._dp
  rec = 2._dp


  ! Read keplerian elements.
  !filename = 'kep_el.txt'
  !call read_kep_el(filename, size(orbits), 1, kep_el)

  ! Evolve keplerian elements per century.
  !call kep_el_teph(kep_el, orbits, 2451545._dp)

  ! Calculate the state vectors.
  !call state_vec(orbits,1d-6)

  ! Evolve the trajectories with time.
  open(unit = 1, file = 'solar_system.dat', status = 'replace', action = 'write', position = 'rewind')
  call integrate_orbits(orbits, start, finish, dt, rec)
  close(1)

  call system('mkdir tmp')
  call pngterm('solar_system', plot_size=[800,800], font='Helvetica', font_size=1)
  !adplot2d(dfilename,pngname,interval,step,using,nplots)
  call adplot3d('solar_system','solar_system',[0,nint(finish/(rec*dt))],1,[(i,i=1,15)],size([(i,i=1,15)])/3)
  call system('gnuplot solar_system.gnu')
  call system('ffmpeg -i tmp/solar_system%d.png animation.mp4')
  call system('rm -rf tmp')


end program solar_system_main
