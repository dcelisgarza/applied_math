module num_diff
  use nrtype
contains
  function sifodi1(func,x,h)
    ! si := simple, fo := forward, di1 = 1st derivative
    implicit none
    real(dp), intent(in) :: x, h ! Variable and step size.
    real(dp)             :: sifodi1
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    div0: if (h == 0) then
      write(*,*) ' Error: sifodi1:: Division by zero.'
      return
    end if div0
    sifodi1 = (func(x+h) - func(x))/h
  end function sifodi1

  function eicedi1(func,x,h)
    ! ei := eigth order, ce := central, di1 := 1st derivative
    implicit none
    real(dp), intent(in) :: x, h ! Variable and step size.
    real(dp)             :: eicedi1
    real(dp)             :: ho8, ho8p6   ! h/8
    real(dp), parameter  :: c4 = 1./280, c3 = 4./105, c2 = 0.2, c1 = 0.8
    interface functn
      function func(x)
        use nrtype
        real(dp), intent(in)           :: x    ! Variables
        real(dp)                       :: func ! Function
      end function func
    end interface functn

    ho8   = 0.125_dp*h
    ho8p6 = ho8**6._dp

    div0: if (ho8p6 == 0) then
      write(*,*) ' Error: eicedi1 :: Division by zero.'
      return
    end if div0
    eicedi1 = (c4*func(x-ho8*4._dp)-c3*func(x-ho8*3._dp)+c2*func(x-ho8*2._dp)&
              -c1*func(x-ho8)&
              +c1*func(x+ho8)&
              -c4*func(x+ho8*4._dp)+c3*func(x+ho8*3._dp)-c2*func(x+ho8*2._dp))/ho8p6
  end function eicedi1

end module num_diff
