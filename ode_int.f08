module ode_int
  implicit none
  contains
  subroutine rk4g(derivs,x,yi,yf,h)
    !================================================================!
    ! Solve n first-order ODEs or n/2 second-order.                  !
    ! Runge-Kutta 4 Gill's method                                    !
    ! Daniel Celis Garza 24 Sept. 2015                               !
    !----------------------------------------------------------------!
    ! f(x,y,dy,n) = ODEs to be solved                                !
    ! n           = ODE system dimension                             !
    ! h           = step size                                        !
    ! x           = independent variable @ step i                    !
    ! yi()        = array of dependent variables @ step i            !
    ! ys()        = array of dependent variables @ stage s of step i !
    ! yf()        = array of dependent variables @ step i+1          !
    ! er()        = array of integration errors                      !
    ! ki()        = array of Runge-Kutta k/h                         !
    ! j           = counter variable                                 !
    !================================================================!
    implicit none
    double precision, intent(in) :: x, h, yi(:)
    double precision, intent(out) :: yf(:)
    double precision, dimension(size(yi)) :: k1, k2, k3, k4, ys
    double precision h2,c1,c2,c3,c4,c5,c6,c7
    parameter (c1=sqrt(2.0),c2=-0.5*c1,c3=2.0-c1,c4=2.0+c1,c5=c2-0.5,c6=0.5*c3,c7=0.5*c4)
    interface
      subroutine derivs(x,y,dydx)
        implicit none
        double precision, intent(in) :: x, y(:)
        double precision, intent(out) :: dydx(:)
      end subroutine
    end interface

    h2=0.5*h
    ! Calculate k1
    call derivs(x, yi, k1) ! Call the derivatives
    ys = yi + h2*k1      ! Go through all equations and prepare for the next stage.

    ! Calculate k2
    call derivs(x+h2, ys, k2)
    ys = yi + h*(c5*k1 + c6*k2)

    ! Calculate k3
    call derivs(x+h2, ys, k3)
    ys = yi + h*(c2*k2 + c7*k3)

    ! Calculate k4 and yf
    call derivs(x+h, ys, k4)
    yf  = yi + h*(k1 + c3*k2 + c4*k3 + k4)/6.0
  end subroutine rk4g

  subroutine rkck(derivs,x,yi,yf,er,h)
    !================================================================!
    ! Solve n first-order ODEs or n/2 second-order.                  !
    ! Cash-Karp RK45                                                 !
    ! Daniel Celis Garza 24 Sept. 2015                               !
    !----------------------------------------------------------------!
    ! f(x,y,dy,n) = ODEs to be solved                                !
    ! n           = ODE system dimension                             !
    ! h           = step size                                        !
    ! x           = independent variable @ step i                    !
    ! yi()        = array of dependent variables @ step i            !
    ! ys()        = array of dependent variables @ stage s of step i !
    ! yf()        = array of dependent variables @ step i+1          !
    ! er()        = array of integration errors                      !
    ! ki()        = array of Runge-Kutta k/h                         !
    ! ci          = Butcher Table c-vector                           !
    ! aij         = Butcher Table A-matrix                           !
    ! bi          = Butcher Table b-vector                           !
    ! dbi         = b-b* vector difference for error calculation     !
    ! j           = counter variable                                 !
    !================================================================!
    implicit none
    double precision, intent(in) :: x, yi(:), h
    double precision, intent(out) :: yf(:), er(:)
    double precision, dimension(size(yi)) :: k1, k2, k3, k4, k5, k6, ys
    double precision c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,&
                     a51,a52,a53,a54,a61,a62,a63,a64,a65,b1,b3,b4,b6,&
                     db1,db3,db4,db5,db6
    parameter ( c2=.2, c3=.3, c4=.6, c5=1., c6=.875, &
                a21=.2, &
                a31=.075, a32=.225, &
                a41=.3, a42=-.9, a43=1.2, &
                a51=-11./54., a52=2.5, a53=-70./27., a54=35./27., &
                a61=1631./55296., a62=175./512., a63=575./13824., a64=44275./110592., a65=253./4096., &
                b1=37./378., b3=250./621., b4=125./594., b6=512./1771., &
                db1=b1-2825./27648., db3=b3-18575./48384., db4=b4-13525./55296., db5=-277./14336., db6=b6-0.25 )
    interface
      subroutine derivs(x,y,dydx)
        implicit none
        double precision, intent(in) :: x, y(:)
        double precision, intent(out) :: dydx(:)
      end subroutine
    end interface

    ! Calculate k1
    call derivs(x, yi, k1)
    ys = yi + h*a21*k1

    call derivs(x+c2*h, ys, k2)
    ys = yi + h*(a31*k1+a32*k2)

    call derivs(x+c3*h, ys, k3)
    ys = yi + h*(a41*k1+a42*k2+a43*k3)

    call derivs(x+c4*h, ys, k4)
    ys = yi + h*(a51*k1+a52*k2+a53*k3+a54*k4)

    call derivs(x+c5*h, ys, k5)
    ys = yi + h*(a61*k1+a62*k2+a63*k3+a64*k4+a65*k5)

    call derivs(x+c6*h, ys, k6)
    yf = yi + h*(b1*k1+b3*k3+b4*k4+b6*k6)

    ! Calculate error
    er = db1*k1+db3*k3+db4*k4+db5*k5+db6*k6
  end subroutine rkck

  subroutine rkcka(derivs,x,y,der,h,hmin)
    !================================================================!
    ! Adaptive-step Cash-Karp RK45                                   !
    ! Daniel Celis Garza 24 Sept. 2015                               !
    !----------------------------------------------------------------!
    ! f(x,y,dy,n) = ODEs to be solved                                !
    ! n           = ODE system dimension                             !
    ! h           = current step size                                !
    ! ht          = trial step size                                  !
    ! hn          = new step size                                    !
    ! dy          = array of derivatives                             !
    ! der         = desired error                                    !
    ! tiny        = prevent division by zero                         !
    ! s           = safety factor for changing h                     !
    ! shrink      = (-1/p-1) where p is the error's order            !
    ! grow        = (-1/p) where p is the error's order              !
    ! ercor       = (5/s)**(1/grow), know when to increase h greatly !
    ! yscal(:)    = scaling the error to the function's range        !
    ! j           = counter variable                                 !
    !================================================================!
    implicit none
    double precision, intent(inout) :: x, y(:), h
    double precision, intent(in) :: der, hmin
    double precision, dimension(size(y)) :: yn, yscal, dydx, er
    double precision local_der, ht, hn, maxer, tiny, s, shrink, grow, ercor
    parameter (tiny=1e-30,s=0.9,shrink=-0.25,grow=-0.2,ercor=1.89e-4)
    interface
      subroutine derivs(x,y,dydx)
        implicit none
        double precision, intent(in) :: x, y(:)
        double precision, intent(out) :: dydx(:)
      end subroutine
    end interface

    local_der = der
    call derivs(x, y, dydx)
    yscal(:) = abs(y(:)) + abs(h*dydx(:)) + tiny
    ! When global errors are a concern use:
    ! yscal(:) = abs(h*dy(:)) + tiny
    do
      call rkck(derivs,x,y,yn,er,h)
      maxer = maxval(abs(er(:)/yscal(:)))/local_der
      if (maxer <= 1.) exit
      ht = s*h*maxer**shrink
      h = sign(max(abs(ht),0.1*abs(h)),h) ! Prevents jumps of more than an order of magnitude.
      if (abs(h)<hmin .or. x+h == x) then
        ! This means that the error is so strict h is effectively zero,
        ! so the program enters an infinite loop.
        ! Increasing h to the minimum value.
        h = hmin
        ! By locally increasing the desired error, we further prevent the infinite loop.
        ! The error goes back to normal for the next integration step.
        local_der = local_der + der
      end if
    end do

    if (maxer > ercor) then
      hn = s*h*maxer**grow
    else
      hn = 5.*h
    end if

    if(abs(hn) < abs(hmin)) hn = hmin

    ! Updating for the next step.
    x = x + h    ! Taking a step.
    h  = hn      ! Setting the next step length.
    y(:) = yn(:) ! Setting all the dependent variables to their new starting points.
  end subroutine rkcka
end module ode_int