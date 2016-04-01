program solar_system_main
  use orbital_mech
  use plot
  implicit none
  real(dp) :: kep_el(20,10) ! Keplerian elements of the bodies of the solar system.
  type(orb_par) :: orbits(10) ! Orbits of the solar system.
  character(:), allocatable :: filename
  integer  :: i, nrec ! counter, number of records
  real(dp) :: t, dt, start, finish, rec! time, increment in time, start, finish, # steps per recording.
  real(dp), allocatable :: dummy(:), ax(:,:), max_ax(:,:), min_ax(:,:)
  real(4) :: maxrange(3), minrange(3)
  real(4) :: xrange(2), yrange(2), zrange(2)
  !type(titles), allocatable :: title(:)
  character(20) :: title(10)
  !character(:), allocatable :: title(:)

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

  ! Calculating the range for the animation of inner solar system.
  nrec = nint(finish/(rec*dt)) + 1
  open(unit = 2, file = 'solar_system.dat', status = 'old', action = 'read', position = 'rewind')
  allocate(ax(15,nrec), max_ax(3,5), min_ax(3,5))
  iranges: do i = 1, nrec
    read(2,*) ax(:,i)
  end do iranges
  close(2)

  min_ax = reshape(minval(ax,dim=2), [3,5])
  max_ax = reshape(maxval(ax,dim=2), [3,5])
  minrange = minval(min_ax,dim=2)
  maxrange = maxval(max_ax,dim=2)
  xrange = [minrange(1), maxrange(1)]
  yrange = [minrange(2), maxrange(2)]
  zrange = [minrange(3), maxrange(3)]
  ! Call gnuplot.
  call system('mkdir tmp')
  call pngterm('isolar_system', plot_size=[1000,1000], font='TeX Gyre Pagella', font_size=13)
  !call range(xrange,yrange,zrange)
  call range([-1.8,1.8],[-1.8,1.8],[-0.06,0.06])
  call plt_labels('x, A.U.', 'y, A.U.', 'z, A.U.', 'Inner Solar System')
  call ticks(xticks=[-1.8,0.4,1.8],mxticks=2,yticks=[-1.8,0.4,1.8],myticks=2,zticks=[-0.06,0.02,0.06], mzticks=2)
  call grid(xticks=[1,0],yticks=[1,0],zticks=[1,0], linestyle=[(0,i=1,6)])
  call xyplane(-0.06)
  title(1) = 'Sun'
  title(2) = 'Mercury'
  title(3) = 'Venus'
  title(4) = 'Earth'
  title(5) = 'Mars'
  call adplot3d('solar_system','isolar_system',[0,nrec],1,[(i,i=1,15)],size([(i,i=1,15)])/3, title)

  ! Calculating the range for the animation of outer solar system.
  nrec = nint(finish/(rec*dt)) + 1
  open(unit = 2, file = 'solar_system.dat', status = 'old', action = 'read', position = 'rewind')
  allocate(dummy(15))
  oranges: do i = 1, nrec
    read(2,*) dummy(:),ax(:,i)
  end do oranges
  close(2)

  min_ax = reshape(minval(ax,dim=2), [3,5])
  max_ax = reshape(maxval(ax,dim=2), [3,5])
  minrange = minval(min_ax,dim=2)
  maxrange = maxval(max_ax,dim=2)
  xrange = [minrange(1), maxrange(1)]
  yrange = [minrange(2), maxrange(2)]
  zrange = [minrange(3), maxrange(3)]
  ! Call gnuplot.
  call pngterm('osolar_system', plot_size=[1000,1000], font='TeX Gyre Pagella', font_size=13)
  call range([-31.,44.0],[-33.,45.],[-15.,9.])
  call plt_labels('x, A.U.', 'y, A.U.', 'z, A.U.', 'Outer Solar System')
  call ticks(xticks=[-31.,15.,44.0],mxticks=2,yticks=[-33.,10.,47.],myticks=2,zticks=[-15.,4.,9.], mzticks=2)
  call grid(xticks=[1,0],yticks=[1,0],zticks=[1,0], linestyle=[(0,i=1,6)])
  call xyplane(-15.)
  title(6) = 'Jupiter'
  title(7) = 'Saturn'
  title(8) = 'Uranus'
  title(9) = 'Neptune'
  title(10) = 'Pluto'
  call adplot3d('solar_system','osolar_system',[0,nrec],5,[(i,i=16,30)],size([(i,i=16,30)])/3, title(6:10))

  call system('gnuplot isolar_system.gnu')
  !call system('gnuplot osolar_system.gnu')
  call system('ffmpeg -i tmp/isolar_system%d.png isolar_system.mp4')
  !call system('ffmpeg -i tmp/osolar_system%d.png osolar_system.mp4')
  call system('rm -rf tmp')
  call system ('make -f make_sol.txt clean')


end program solar_system_main
