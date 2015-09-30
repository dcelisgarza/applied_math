module pend_len
! Allows the declaration of global variables whose values can be assigned in other places.
double precision :: l
end module

module dpend
double precision :: l1, l2, m1, m2, ms
end module

module gonewiththewind
double precision :: as, ar, k
end module

module per_harv
double precision :: a, h, r
end module

program main
use pend_len
use dpend
use gonewiththewind
use per_harv
implicit none
! call simple_pendulum
 call double_pendulum
 call double_pendulum_adapt
! call gone_with_the_wind
! call periodic_harvesting
! call periodic_harvesting2
! call periodic_harvesting_adapt
end program main

subroutine periodic_harvesting
! Excercise 7.2.19 of non-linear dynamics and chaos Strogatz
use per_harv
implicit none
integer, parameter :: n = 1
double precision ti, tf, dt, tmax
double precision xi(n), xf(n)
integer i
external pharv

open(unit=1, file='pharv.dat')
! initial data
ti    = 0.0
tmax  = 100.0
dt    = 0.01
a     = 0.5
r     = 1.0
h     = 0.2  ! criticals are roughly h<r/4, h<r/(4*(1+a))
xi(1) = 1.0  ! dimensionless model, the maximum population is really a proportion.

! print initial conditions
write(1,*) ti, xi(1)

! integration of ODEs
do while(ti <= tmax)
   tf = ti + dt
   call rk4n(pharv,ti,tf,xi,xf,n)
   write(1,*) tf, xf(1)

! preparing for next step
   ti = tf
   do i = 1,n
      xi(i) = xf(i)
   end do
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
end subroutine periodic_harvesting

subroutine periodic_harvesting2
! Excercise 7.2.19 of non-linear dynamics and chaos Strogatz
use per_harv
use functions
implicit none
integer, parameter :: n = 1
double precision ti, tf, dt, tmax
double precision xi(n), xf(n), er(n)
integer i
external pharv

open(unit=1, file='pharv2.dat')
! initial data
ti    = 0.0
tmax  = 100.0
dt    = 0.01
a     = 0.5
r     = 1.0
h     = 0.2  ! criticals are roughly h<r/4, h<r/(4*(1+a))
xi(1) = 1.0  ! dimensionless model, the maximum population is really a proportion.

! print initial conditions
write(1,*) ti, xi(1)

! integration of ODEs
do while(ti <= tmax)
   !call rkckn(pharv,ti,xi,xf,er,dt,n)
   tf = ti + dt
   write(1,*) tf, xf(1)

! preparing for next step
   ti = tf
   do i = 1,n
      xi(i) = xf(i)
   end do
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
end subroutine periodic_harvesting2

subroutine periodic_harvesting_adapt
! Excercise 7.2.19 of non-linear dynamics and chaos Strogatz
use per_harv
use functions
implicit none
integer, parameter :: n = 1
double precision ti, tf, dt, tmax
double precision xi(n), xf(n), er(n)
double precision der, dtmin
integer i
external pharv

open(unit=1, file='pharv_adapt.dat')
! initial data
ti    = 0.0
tmax  = 100.0
dt    = 0.01
der   = 1e-7
dtmin = 1e-4
a     = 0.5
r     = 1.0
h     = 0.2  ! criticals are roughly h<r/4, h<r/(4*(1+a))
xi(1) = 1.0  ! dimensionless model, the maximum population is really a proportion.

! print initial conditions
write(1,*) ti, xi(1)

! integration of ODEs
do while(ti <= tmax)
  !call rkcka(pharv,ti,tf,xi,xf,der,dt,dtmin,n)
  write(1,*) tf, xf(1)
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
end subroutine periodic_harvesting_adapt

subroutine pharv(t,x,dx,n)
! Excercise 8.5.4 strogatz
use per_harv
implicit none
integer n
double precision t, x(n), dx(n)

dx(1) = r*x(1)*(1-x(1)) - h*(1 + a*sin(t))
end subroutine pharv

subroutine gone_with_the_wind
! Excercise 7.2.19 of non-linear dynamics and chaos Strogatz
use gonewiththewind
implicit none
integer, parameter :: n = 2
double precision ti, tf, dt, tmax
double precision xi(n), xf(n)
integer i
external gwtw

open(unit=1, file='gwtw.dat')
! initial data
ti    = 0.0
tmax  = 100.0
dt    = 0.01
as    = 1.2
ar    = 1.0
k     = 15.0
xi(1) = 0.0
xi(2) = 0.0

! print initial conditions
write(1,*) ti, xi(1), xi(2)

! integration of ODEs
do while(ti <= tmax)
   tf = ti + dt
   call rk4n(gwtw,ti,tf,xi,xf,n)
   write(1,*) tf, xf(1), xf(2)

! preparing for next step
   ti = tf
   do i = 1,n
      xi(i) = xf(i)
   end do
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
end subroutine gone_with_the_wind

subroutine gwtw(t,x,dx,n)
!======================================================!
! Daniel Celis Garza  13 Dec. 2014                     !
! Love affair between Scarlett O'Hara and Rhett Butler !
! x(1) = r    dx(1) = r'                               !
! x(2) = s    dx(2) = s'                               !
!======================================================!
use gonewiththewind
implicit none
integer n
double precision t, x(n), dx(n)

dx(1) = -x(1) + as + k*x(2)*exp(-x(2))
dx(2) = -x(2) + ar + k*x(1)*exp(-x(1))
end subroutine gwtw

subroutine double_pendulum
use dpend
use functions
implicit none
integer, parameter :: n = 4
double precision, parameter :: rad = 3.1415926/180.0
double precision ti, tf, dt, tmax, thsi
double precision xi(n), xf(n), dxi(n),er(n)
double precision start, finish
integer i
interface
  subroutine dpen(x,y,dx)
    implicit none
    double precision, intent(in) :: x, y(:)
    double precision, intent(out) :: dx(:)
  end subroutine
end interface

open(unit=1, file = 'd_pendulum.dat')
! initial data
ti     = 0.0
tmax   = 50.0
dt     = 5.0e-3
l1     = 5.0
l2     = 9.0
m1     = 13.0
m2     = 50.0
xi(1)  = 173.0*rad
xi(2)  = -86.0*rad
dxi(1) = -3.0*rad
dxi(2) = 35.0*rad
ms     = m1 + m2
thsi   = xi(1) - xi(2)
xi(3)  = ms*l1*l1*dxi(1) + m2*l1*l2*dxi(2)*cos(thsi)
xi(4)  = m2*l2*l2*dxi(2) + m2*l1*l2*dxi(1)*cos(thsi)

! print initial conditions
write(1,*) ti,xi(1),xi(2),xi(3),xi(4),l1*sin(xi(1)),-l1*cos(xi(1)),l1*sin(xi(1))+l2*sin(xi(2)),-l1*cos(xi(1))-l2*cos(xi(2))
CALL CPU_TIME(start)
! integration of ODEs
do while(ti <= tmax)
   tf = ti + dt
   call rkckn(dpen,ti,xi,xf,er,dt)
   !call rkckn(dpen,ti,xi,xf,er,dt,n)
   write(1,*) tf,xf(1),xf(2),xf(3),xf(4),l1*sin(xf(1)),-l1*cos(xf(1)),l1*sin(xf(1))+l2*sin(xf(2)),-l1*cos(xf(1))-l2*cos(xf(2))

! prepare for next step
   ti = tf
   do i = 1,n
      xi(i) = xf(i)
      end do
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
CALL CPU_TIME(finish)
print '("Time = ",f6.3," seconds.")',finish-start
end subroutine double_pendulum

subroutine double_pendulum_adapt
use dpend
use functions
implicit none
integer, parameter :: n = 4
double precision, parameter :: rad = 3.1415926/180.0
double precision ti, tf, dt, dtmin, der, tmax, thsi, t
double precision x(n)
double precision start, finish
integer i
interface
  subroutine dpen(x,y,dx)
    implicit none
    double precision, intent(in) :: x, y(:)
    double precision, intent(out) :: dx(:)
  end subroutine
end interface
double precision, dimension(size(x)) :: xi, xf, dxi, dx

open(unit=1, file = 'd_pendulum_adapt.dat')
! initial data
ti     = 0.0
tmax   = 50.0
dt     = 5.0e-3
dtmin  = 4.375e-3
der    = 1.3125e-7
l1     = 5.0
l2     = 9.0
m1     = 13.0
m2     = 50.0
xi(1)  = 173.0*rad
xi(2)  = -86.0*rad
dxi(1) = -3.0*rad
dxi(2) = 35.0*rad
ms     = m1 + m2
thsi   = xi(1) - xi(2)
xi(3)  = ms*l1*l1*dxi(1) + m2*l1*l2*dxi(2)*cos(thsi)
xi(4)  = m2*l2*l2*dxi(2) + m2*l1*l2*dxi(1)*cos(thsi)
! print initial conditions
write(1,*) ti,xi(1),xi(2),xi(3),xi(4),l1*sin(xi(1)),-l1*cos(xi(1)),l1*sin(xi(1))+l2*sin(xi(2)),-l1*cos(xi(1))-l2*cos(xi(2))
t = ti
x(:) = xi(:)
CALL CPU_TIME(start)
! integration of ODEs
do while(t <= tmax)
  call rkcka(dpen,t,x,der,dt,dtmin)
  !call rkcka(dpen,ti,tf,xi,xf,der,dt,dtmin,n)
  write(1,*) t,x(1),x(2),x(3),x(4),l1*sin(x(1)),-l1*cos(x(1)),l1*sin(x(1))+l2*sin(x(2)),-l1*cos(x(1))-l2*cos(x(2))
end do
CALL CPU_TIME(finish)
print '("Time = ",f6.3," seconds.")',finish-start
end subroutine double_pendulum_adapt

subroutine dpen(t,x,dx)
! x(1) = theta1    dx(1) = theta1'
! x(2) = theta2    dx(2) = theta2'
! x(3) = p1        dx(3) = p1'
! x(4) = p2        dx(4) = p2'
use dpend
implicit none
double precision t, ths, den, c1, c2
double precision, intent(inout) :: x(:)
double precision, intent(out):: dx(:)
double precision, parameter :: g = 9.81
ths   = x(1) - x(2)
den   = m1 + m2*sin(ths)*sin(ths)
c1    = (x(3)*x(4)*sin(ths))/(l1*l2*den)
c2    = (l2*l2*m2*x(3)*x(3) + l1*l1*ms*x(4)*x(4) - l1*l2*m2*x(3)*x(4)*cos(ths))/(2*l1*l1*l2*l2*den*den)*sin(2*ths)
dx(1) = (l2*x(3)-l1*x(4)*cos(ths))/(l1*l1*l2*den)
dx(2) = (l1*ms*x(4)-l2*m2*x(3)*cos(ths))/(l1*l2*l2*m2*den)
dx(3) = -ms*g*l1*sin(x(1)) - c1 + c2
dx(4) = -m2*g*l2*sin(x(2)) + c1 - c2
end subroutine dpen

subroutine simple_pendulum
!=======================================================!
! Simple pendulum routine                               !
!=======================================================!
use pend_len  ! Places the module into the subroutine along the variable values as assigned by the module of main program.
implicit none
integer, parameter :: n = 2
double precision, parameter :: rad = 3.1415926/180.0
double precision ti, tf, dt, tmax
double precision xi(n), xf(n)
integer i
external spen   ! external allows other functions to be used as arguments.

open(unit=1, file='pendulum.dat')
! initial data
ti    = 0.0     ! initial time
tmax  = 1.0     ! maximum time
dt    = 0.01    ! step size
l     = 1.0     ! assigns the global value for l, it applies to any subroutines or functions called by this subroutine.
xi(1) = 90*rad  ! initial angle
xi(2) = 0.0*rad ! initial angular speed

! print initial conditions
write(1,102) ti, xi(1), xi(2), l*sin(xi(1)), -l*cos(xi(1))

! integration of ODEs
do while(ti<=tmax)
   tf = ti + dt
   call rk4n(spen,ti,tf,xi,xf,n)
   write(1,102) tf, xf(1), xf(2), l*sin(xf(1)), -l*cos(xf(1))
! prepare next step
   ti = tf
   do i = 1, n
      xi(i) = xf(i)
   end do
end do

100 format(5x,'t',11x,'x',11x,'y',11x,'dx/dt',7x,'dy/dt')
102 format(5(1pe12.3))
end subroutine simple_pendulum

subroutine spen(t,x,dx,n)
use pend_len
!=======================================================!
! System of ODEs for the rk4n subroutine.
! Daniel Celis Garza 11 Dec. 2014
!-------------------------------------------------------!
! Simple pendulum
! x(1) = theta         dx(1) = theta'
! x(2) = theta'        dx(2) = theta''
!=======================================================!
implicit none
integer n
double precision t, x(n), dx(n)
double precision, parameter :: g = 9.81                 ! acceleration in m/s^2 and string length
! first order
dx(1) = x(2)
! second order
dx(2) = (-g/l) * sin(x(1))
end subroutine spen

subroutine rk4n(f,ti,tf,xi,xf,n)
!=======================================================!
! Solve n first-order ODEs or n/2 second-order.         !
! Runge-Kutta 4                                         !
! Daniel Celis Garza 10 Dec. 2014                       !
!-------------------------------------------------------!
! f(t,x,dx,n) = parametrised equations of motion        !
! ti          = independent variable @ step i           !
! tf          = independent variable @ step i+1         !
! xi()        = array of dependent variables @ step i   !
! xf()        = array of dependent variables @ step i+1 !
! n           = # of parametrised equations of motion   !
!-------------------------------------------------------!
! j           = counter variable                        !
! h           = time step                               !
! t           = current time                            !
! ki()        = array of Runge-Kutta k's                !
!=======================================================!
implicit none
integer n, j
double precision ti, tf, h, t, xi(n), xf(n), k1(n), k2(n), k3(n), k4(n), x(n), dx(n)

h = tf - ti ! This will be replaced by the error-correcting routine.
t = ti

! Calculate k1
call f(t,xi,dx,n)            ! Call the equations of motion
do j = 1, n                  ! Go through all equations of motion
   k1(j) = h*dx(j)           ! Calculate the value of k1 for each equation of motion
   x(j)  = xi(j) + k1(j)/2.0 ! Calculate the next value of x for each equation (to be used in the next ki)
end do

! Calculate k2
call f(t+h/2.0,x,dx,n)
do j = 1, n
   k2(j) = h*dx(j)
   x(j)  = xi(j) + k2(j)/2.0
end do

! Calculate k3
call f(t+h/2.0,x,dx,n)
do j = 1, n
   k3(j) = h*dx(j)
   x(j)  = xi(j) + k3(j)/2.0
end do

! Calculate k4 and xf
call f(t+h,x,dx,n)
do j = 1, n
   k4(j) = h*dx(j)
   xf(j)  = xi(j) + k1(j)/6.0 + k2(j)/3.0 + k3(j)/3.0 + k4(j)/6.0
end do
end subroutine rk4n

!subroutine rk4n_adapt(f,ti,tf,xi,xf,n)
!=======================================================
! Dynamic adjustment of h.
! Daniel Celis 10 Dec. 2014
! ert = 0 ; absolute error
! ert = 1 ; relative error
! ert = 2 ; mixed error
!=======================================================
!integer j, n, ert
!double precision ti, tf, h, t, xi(n), xf(n), k1(n), k2(n), k3(n), k4(n), x(n), dx(n)

!end subroutine rk4nvh
