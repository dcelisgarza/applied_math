program attractors
  use ode
  use dyn_sys
  implicit none
  real(dp) :: t, start, finish, step, min_step, desired_error
  real(dp) :: y(3)
  integer  :: i, nstepwrite

  open(unit = 1, file = "parameters.dat")

  ! Chen
  ! Read integration parameters
  read(1,*) start, finish, step, min_step, desired_error, nstepwrite
  read(1,*) a, c, b
  read(1,*) y
  open(unit = 2, file = "chen.dat")

  ! Write initial conditions.
  write(2,*) y
  t = 0._dp
  i = 0

  ! Integrate
  dchen: do while (t < finish)
    call rkcka(chen, t, y, desired_error, step, min_step)
    i = i + 1
    if(mod(i, nstepwrite) == 0) write(2,*) y
  end do dchen
  close(2)

  ! Lu Chen
  ! Read integration parameters.
  read(1,*) start, finish, step, min_step, desired_error, nstepwrite
  read(1,*) a, c, b, u
  read(1,*) y
  open(2, file = "lu_chen.dat")

  write(2,*) y
  t = 0._dp
  i = 0
  dlu_chen:do while (t < finish)
    call rkcka(lu_chen, t, y, desired_error, step, min_step)
    i = i + 1
    if(mod(i, nstepwrite) == 0) write(2,*) y
  end do dlu_chen
  close(2)


  ! Modified Chaotic Chua
  ! Read integration parameters.
  read(1,*) start, finish, step, min_step, desired_error, nstepwrite
  read(1,*) alpha, beta, a, b, c, d
  read(1,*) y
  open(2, file = "mc_chua.dat")

  write(2,*) y
  t = 0._dp
  i = 0
  dmc_chua:do while (t < finish)
    call rkcka(mc_Chua, t, y, desired_error, step, min_step)
    i = i + 1
    if(mod(i, nstepwrite) == 0) write(2,*) y
  end do dmc_chua
  close(2)


end program attractors
