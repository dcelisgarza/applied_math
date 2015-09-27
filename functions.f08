module functions
implicit none
contains
  function heaviside(x)
    !=======================================================!
    ! Heaviside function                                    !
    ! Daniel Celis Garza 15 Jan 2015                        !
    !=======================================================!
    double precision heaviside, x
    heaviside = 0.0
    if (x .ge. 0.0) heaviside = 1.0
  end function heaviside

  function kdelta(i,j)
    !=======================================================!
    ! Kroencker delta function                              !
    ! Daniel Celis Garza 15 Jan 2015                        !
    ! delta_ij                                              !
    !=======================================================!
    integer i, j
    double precision kdelta
    kdelta = 0.0
    if (i .eq. j) kdelta = 1.0
  end function kdelta

  subroutine rk4gn(f,x,yi,yf,h,n)
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
    integer n, j
    double precision x,dy(n),yi(n),ys(n),yf(n),&
                     k1(n),k2(n),k3(n),k4(n),&
                     h,h2,c1,c2,c3,c4,c5,c6,c7
    parameter (c1=sqrt(2.0),c2=-0.5*c1,c3=2.0-c1,c4=2.0+c1,c5=c2-0.5,c6=0.5*c3,c7=0.5*c4)

    h2=0.5*h
    ! Calculate k1
    call f(x, yi, k1, n)          ! Call the derivatives
    do j = 1, n                   ! Go through all equations of motion
      ys(j) = yi(j) + h2*k1(j)    ! Calculate the next value of x for each equation
    end do                        ! (to be used in the next ki)

    ! Calculate k2
    call f(x+h2, ys, k2, n)
    do j = 1, n
       ys(j) = yi(j) + h*(c5*k1(j) + c6*k2(j))
    end do

    ! Calculate k3
    call f(x+h2, ys, k3, n)
    do j = 1, n
       ys(j) = yi(j) + h*(c2*k2(j) + c7*k3(j))
    end do

    ! Calculate k4 and yf
    call f(x+h, ys, k4, n)
    do j = 1, n
       yf(j)  = yi(j) + h*(k1(j) + c3*k2(j) + c4*k3(j) + k4(j))/6.0
    end do
  end subroutine rk4gn

  subroutine rkckn(f,x,yi,yf,er,h,n)
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
    integer n, j
    double precision x,dy(n),yi(n),ys(n),yf(n),&
                     k1(n),k2(n),k3(n),k4(n),k5(n),k6(n),er(n),&
                     h,c2,c3,c4,c5,c6,a21,a31,a32,a41,a42,a43,&
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

    ! Calculate k1
    call f(x, yi, k1, n)
    do j = 1, n
      ys(j) = yi(j) + h*a21*k1(j)
    end do

    call f(x+c2*h, ys, k2, n)
    do j = 1, n
      ys(j) = yi(j) + h*(a31*k1(j)+a32*k2(j))
    end do

    call f(x+c3*h, ys, k3, n)
    do j = 1, n
      ys(j) = yi(j) + h*(a41*k1(j)+a42*k2(j)+a43*k3(j))
    end do

    call f(x+c4*h, ys, k4, n)
    do j = 1, n
      ys(j) = yi(j) + h*(a51*k1(j)+a52*k2(j)+a53*k3(j)+a54*k4(j))
    end do

    call f(x+c5*h, ys, k5, n)
    do j = 1, n
      ys(j) = yi(j) + h*(a61*k1(j)+a62*k2(j)+a63*k3(j)+a64*k4(j)+a65*k5(j))
    end do

    call f(x+c6*h, ys, k6, n)
    do j = 1, n
      yf(j) = yi(j) + h*(b1*k1(j)+b3*k3(j)+b4*k4(j)+b6*k6(j))
    end do

    ! Calculate error
    do j = 1, n
      er(j) = db1*k1(j)+db3*k3(j)+db4*k4(j)+db5*k5(j)+db6*k6(j)
    end do
  end subroutine rkckn

  subroutine rkcka(f,xi,xf,yi,yf,der,h,hmin,n)
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
    ! yscal(j)    = scaling the error to the function's range        !
    ! j           = counter variable                                 !
    !================================================================!
    external f
    integer n, j
    double precision xi,xf,xtest,yi(n),yf(n),dy(n),yscal(n),der,er(n),h,ht,hn,hmin,maxer,miner,&
                     tiny,s,shrink,grow,ercor
    parameter (tiny=1e-30,s=0.9,shrink=-0.25,grow=-0.2,ercor=1.89e-4)
    call f(xi, yi, dy, n)
    do j = 1, n
       yscal(j) = abs(yi(j)) + abs(h*dy(j)) + tiny
      ! When global errors are a concern use:
      ! yscal(j) = abs(h*dy(j)) + tiny
    end do

1   call rkckn(f,xi,yi,yf,er,h,n)
    maxer = 0.
    do j = 1, n
      maxer = max(maxer,abs(er(j)/yscal(j)))
    end do
    maxer = maxer/der

    if (maxer > 1.) then
      ht = s*h*maxer**shrink
      h = max(ht,0.1*h) ! Prevents jumps of more than an order of magnitude.
      xtest = xi+h
      if (xtest == xi) then
        ! This means that the error is so strict h is effectively zero,
        ! so the program enters an infinite loop.
        ! Increasing h to the minimum value,
        h = hmin
        ! and locally increasing the desired error, we prevent this from occurring.
        ! The error goes back to normal for the next integration step.
        der = der + .0001*der
      endif
      go to 1
    else
      if (maxer > ercor) then
        hn = s*h*maxer**grow
      else
        hn = 5.*h
      end if
    endif

    ! Updating for the next step.
    xf = xi + h   ! Taking a step.
    h  = hn       ! Setting the next step length.
    xi = xf       ! Setting the dependent variable as the new starting point.
    do j = 1, n
      yi(j) = yf(j) ! Setting all the dependent variables to their new starting points.
    end do

    ! Making sure we don't go under the minimum value of h.
    !if (h < hmin) then
      ! This way of dealing with it is one of many, hmin could be made smaller automatically.
      ! For example:
      ! hmin = h
      !write(*,*) 'Step size smaller than the minimum specified (', hmin, ').'
      !write(*,*) 'Please define a new minimum value.'
      !write(*,*) 'H_min = '
      !read *, hmin
    !endif
  end subroutine rkcka
end module functions
