program solar_system_main
  use orbital_mech
  implicit none
  type(orb_par) :: planets(1)
  real(dp) :: dt, mu, bom, som, i, ma, eu(3), reu(6)
  integer  :: j

  dt = 8640._dp
  mu = 1.32712440018*10.**20._dp!/149597870700._dp

  ! orb_par(a, ea, ec, ma, mu, ta, dt, n, r, o(4), x(6), eu(3), reu(6))
  eu(1) = 48.33076593_dp*pi/180._dp ! Big omega
  eu(2) = 7.00497902_dp*pi/180._dp  ! I
  eu(3) = (77.45779628_dp-48.33076593_dp)*pi/180._dp ! Small omega

  reu(1) = cos(eu(1))*cos(eu(3)) - sin(eu(1))*cos(eu(2))*sin(eu(3))
  reu(2) = sin(eu(1))*cos(eu(3)) + cos(eu(1))*cos(eu(2))*sin(eu(3))
  reu(3) = sin(eu(2))*sin(eu(3))
  reu(4) = -cos(eu(1))*sin(eu(3)) - sin(eu(1))*cos(eu(2))*cos(eu(3))
  reu(5) = -sin(eu(1))*sin(eu(3)) + cos(eu(1))*cos(eu(2))*cos(eu(3))
  reu(6) = sin(eu(2))*sin(eu(1))
  ma  = (252.2503235_dp-77.45779628_dp)*pi/180._dp

  planets(1) = orb_par(0.38709927_dp, 0._dp, 0.20563593_dp, &
                       ma, mu, &
                       0._dp, 0._dp, sqrt(mu/(0.38709927_dp*0.38709927_dp*0.38709927_dp))&
                       , 0._dp, [0._dp,0._dp,0._dp], [0._dp,0._dp,0._dp,0._dp,0._dp,0._dp]&
                       ,&
                       eu, reu)
 open(unit=1,file='mercury.dat')
 call eorb(planets,10**(-6._dp))

 write(1,*) planets(1)%x(1), planets(1)%x(2)
 planets(1) % dt = dt

 do j = 1, 20
   call eorb(planets,10**(-6._dp))
    write(1,*) planets(1)%x(1), planets(1)%x(2)
 end do
 close(1)
end program solar_system_main
