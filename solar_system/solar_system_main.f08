program solar_system_main
  use orbital_mech
  use plot
  implicit none
  real(dp) :: kep_el(20,10) ! Keplerian elements of the bodies of the solar system.
  type(orb_par) :: orbits(10) ! Orbits of the solar system.
  character(:), allocatable :: filename
  integer  :: i, nrec ! counter, number of records
  real(dp) :: t, dt, start, finish, rec! time, increment in time, start, finish, # steps per recording.
  real(dp), allocatable :: ax(:,:)
  real(dp):: dummy(12)
  real(4) :: maxrange(3), minrange(3)
  real(4) :: xrange(2), yrange(2), zrange(2)
  !type(titles), allocatable :: title(:)
  !character(20) :: title(5)
  character(:), allocatable :: title(:)

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

  ! Calculating the range for the animation.
  nrec = nint(finish/(rec*dt)) + 1
  open(unit = 2, file = 'solar_system.dat', status = 'old', action = 'read', position = 'rewind')
  allocate(ax(3,nrec))
  ranges: do i = 1, nrec
    read(2,*) dummy, ax(:,i)
  !  read(2,*) dummy, ax
  end do ranges
  minrange = minval(ax,dim=2)
  maxrange = maxval(ax,dim=2)
  xrange = [minrange(1), maxrange(1)]
  yrange = [minrange(2), maxrange(2)]
  zrange = [minrange(3)-0.01, maxrange(3)]
  close(2)
  ! Call gnuplot.
  call system('mkdir tmp')
  call pngterm('solar_system', plot_size=[800,800], font='Helvetica', font_size=12)
  call range(xrange,yrange,zrange)
  call format('X, A.U.', 'Y, A.U.', 'Z, A.U.', 'Inner Solar System')
  allocate( character(20) :: title(5) )
!  allocate( character :: title(5)%title )


  !title(1:5) % length = [len('Sun'), len('Mercury'), len('Venus'), len('Earth'), len('Mars')]
  title(1) = 'Sun'
  title(2) = 'Mercury'
  title(3) = 'Venus'
  title(4) = 'Earth'
  title(5) = 'Mars'
  !title(1) % title(1: title(1) % length) = 'Sun'
  !title(2) % title(1: title(2) % length) = 'Mercury'
  !title(3) % title(1: title(3) % length) = 'Venus'
  !title(4) % title(1: title(4) % length) = 'Earth'
  !title(5) % title(1: title(5) % length) = 'Mars'
  !print*, title%leng
  !print*, title(1) % title(1: title(1) % length)


  call adplot3d('solar_system','solar_system',[0,nrec],1,[(i,i=1,15)],size([(i,i=1,15)])/3, title)
  call system('gnuplot solar_system.gnu')
  call system('ffmpeg -i tmp/solar_system%d.png animation.mp4')
  call system('rm -rf tmp')


end program solar_system_main
