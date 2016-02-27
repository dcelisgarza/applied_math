program solar_system_main
  use orbital_mech
  implicit none
  real(dp) :: kep_el(12,9)
  character(:), allocatable :: filename
  integer  :: i
  filename = 'kep_el.txt'
  call read_kep_el(filename, 9, kep_el)
  do i = 1, 9
    write(*,*) kep_el(1:6,i)
    write(*,*) kep_el(7:12,i)
  end do
end program solar_system_main
