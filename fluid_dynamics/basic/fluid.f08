module fluid
  use interpolation, only : cctrilinint
  use fluid_types
  implicit none

  ! Based on
  ! http://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid_flow_for_the_rest_of_us.pdf
contains
  subroutine timestep(kcfl, cell_dim, max_vel, h)
    ! Calculate the time step.
    implicit none
    real(dp), intent(in)    :: kcfl, cell_dim  ! Scaling parameter, cell dimension.
    real(dp), intent(in)    :: max_vel ! Maximum velocity in the grid.
    real(dp), intent(inout) :: h     ! Time step

    ! Time step should be small enough that the maximum velocity field is less than the cell width.
    h = kcfl * cell_dim / max_vel
  end subroutine timestep

  subroutine find_vmax(grid, max_vel)
    implicit none
    type(FluidCell), intent(in) :: grid(:,:,:)
    real(dp), intent(out)       :: max_vel ! Max velocity.
    real(dp) :: tmp_vel ! Trial velocity.
    integer  :: i, j, k

    max_vel = 0.

    dok: do k = 1, size(grid, dim = 3)
      doj: do j = 1, size(grid, dim = 2)
        doi: do i = 1, size(grid, dim = 1)
          ! If the cell type does not contain liquid, skip to the next i.
          if ( grid(i, j, k) % marker == 0 ) cycle doi
          tmp_vel = norm2(grid(i, j, k) % velocity)
          if (tmp_vel > max_vel) max_vel = tmp_vel
        end do doi
      end do doj
    end do dok
  end subroutine find_vmax

  subroutine update_grid(kcfl, grid, grid_size)
    implicit none
    real(dp), intent(in)           :: kcfl
    integer, intent(in)            :: grid_size(:) ! Sizes of the grid. grid_size(1) = size_x, grid_size(2) = size_y, grid_size(3) = size_z
    type(FluidCell), intent(inout) :: grid(:,:,:)
    integer :: i, j, k, l, two_step_counter ! Counters
    integer :: buffer_idx, cell_idx(6) ! Index for buffer at i,j,k, index for the cell at which a buffer will be created.

    ! Set layer of all cells including buffer cells to -1
    grid % layer = -1
    do i = 1, 6
      grid % buffer(i) % layer = -1
    end do

    ! TO DO: MAKE THE PARTICLE MARKERS A LINKED LIST WITH THE INDICES OF THE CELL IN WHICH THEY ARE FOUND. THIS Tvector in each grid cell. We found this be easier and more memory eRIPLE LOOP SHIT CAN GO TO HELL, AND YOU CAN SIMPLY TRAVERSE THE LINKED LIST AND OBTAINING THE INDICES OF THE CORRESPONDING GRID CELL. THIS CAN HELP WHEN EXPANDING THE CODE TO INCLUDE SINKS AND SOURCES.

    ! Update cells that currently have fluid in them.
    dok: do k = 1, grid_size(3)
      doj: do j = 1, grid_size(2)
        doi: do i = 1, grid_size(1)
          ! If the cell type does not contain liquid, skip to the next i.
          if ( grid(i,j,k) % marker == 0 ) cycle doi
          ! If Cell Type is Not a Solid.
          ctns: if ( grid(i,j,k) % cell_type /= 3 ) then
            ! Set cell type to fluid.
            grid(i,j,k) % cell_type = 2
            ! Set layer to zero
            grid(i,j,k) % layer = 0
          end if ctns
        end do doi
      end do doj
    end do dok
    ! Create buffer zone.
    dol: do l = 1, max(2, ceiling(kcfl))

      dok2: do k = 1, grid_size(3)
        doj2: do j = 1, grid_size(2)
          doi2: do i = 1, grid_size(1)

            ! If the grid does NOT contain liquid or air, and its layer is not i-1, cycle.
            if (grid(i,j,k) % cell_type == 3 .or. grid(i,j,k) % layer /= l-1) cycle doi2

            ! Set buffer index to 0.
            buffer_idx = 0
            ! Find relevant cell indices.
            cell_idx  = [i-1, i+1, j-1, j+1, k-1, k+1]
            ! Set two step counter to 0.
            two_step_counter = 0

            ! For each of the six neighbours of grid(i,j,k), update buffer.
            dob: do buffer_idx = 1, 6
              ! Increase two_step_counter every odd step in buffer_idx.
              if ( mod(buffer_idx,2) /= 0 ) two_step_counter = two_step_counter + 1
              ! Increase the buffer index.

              ! Check Buffer is Not Initial.
              bni: if (grid(i,j,k) % buffer(buffer_idx) % initial .eqv. .false.) then

                ! Check if buffer's layer is not -1, and that it is not a solid.
                blbt: if( grid(i,j,k) % buffer(buffer_idx) % layer     == -1 .and. &
                grid(i,j,k) % buffer(buffer_idx) % cell_type /= 3         ) then
                ! Set cell type of buffer to air.
                grid(i,j,k) % buffer(buffer_idx) % cell_type = 1
                grid(i,j,k) % buffer(buffer_idx) % layer     = l
              end if blbt
            else bni
              grid(i,j,k) % buffer(buffer_idx) % initial   = .false.
              grid(i,j,k) % buffer(buffer_idx) % layer     = l
              ! Check whether the cell index which corresponds to the buffer is Out of Bounds.
              ob: if (cell_idx(buffer_idx) < 1 .or. cell_idx(buffer_idx) > grid_size(two_step_counter)) then
                grid(i,j,k) % buffer(buffer_idx) % cell_type = 3
              else ob
                grid(i,j,k) % buffer(buffer_idx) % cell_type = 1
              end if ob
            end if bni
          end do dob

        end do doi2
      end do doj2
    end do dok2

  end do dol

  dok3: do k = 1, grid_size(3)
    doj3: do j = 1, grid_size(2)
      doi3: do i = 1, grid_size(1)
        dom2: do buffer_idx = 1, 6
          if (grid(i,j,k) % buffer(buffer_idx) % layer == -1) grid(i,j,k) % buffer(buffer_idx) % initial = .true.
        end do dom2
      end do doi3
    end do doj3
  end do dok3

end subroutine update_grid

subroutine advect(grid)
  implicit none
  type(FluidCell), intent(in) :: grid(:,:,:)

  !call rkck(derivs, x, yi, yf, er, h)
end subroutine advect

subroutine evolve_velocity(grid, xi, t, h, cell_volume, velocity_offset)
  implicit none
  type(FluidCell), intent(inout) :: grid(:,:,:)
  real(dp), intent(in)           :: xi(:), t, h, cell_volume, velocity_offset(:,:)
  integer :: idx(size(xi))

  ! Set up the index, to save execution time.
  idx = floor(xi)




end subroutine evolve_velocity

subroutine interpolate_velocity(grid, x, idx, cell_volume, velocity_offset)
  ! Interpolate the velocity within the grid from the 8 neighbours.
  implicit none
  type(FluidCell), intent(inout) :: grid(:,:,:) ! Grid
  real(dp), intent(in)           :: x(:), cell_volume, velocity_offset(:,:) ! Coordinates, cell volume
  integer, intent(in)            :: idx(:)
  real(dp) :: tmp_x(size(x))
  integer  :: i ! counter

  interp_vel: do i = 1, size(x)
    ! Correct position offset.
    tmp_x = x - velocity_offset(i,:)

    ! Interpolate velocity.
    grid(idx(1),idx(2),idx(3)) % tmp_velocity(i) = cctrilinint(grid % tmp_velocity(i), tmp_x, idx, cell_volume)
  end do interp_vel

end subroutine interpolate_velocity

subroutine rkck(interp_vel, grid, xi, idx, t, h, cell_volume, velocity_offset)
  !=========================================================!
  ! Backward particle trace                                 !
  ! Cash-Karp RK45                                          !
  ! Daniel Celis Garza 6 June 2016                          !
  !---------------------------------------------------------!
  ! f(x,y,dydx) = ODEs to be solved                         !
  ! Let n := size(y) then for a k'th order system           !
  ! y(1+i*n/k:(i+1)*n/k) := i'th order derivative           !
  !---------------------------------------------------------!
  ! Inputs:                                                 !
  ! derivs = derivatives                                    !
  ! t      = time @ step i                                  !
  ! h      = step size                                      !
  ! cell_volume = cell_volume                               !
  ! velocity_offset(:,:) = offset for the velocities        !
  ! x(:) = position of the particle.                        !
  !---------------------------------------------------------!
  ! Input-Output                                            !
  ! grid(:,:,:) = Fluid dynamics grid                       !
  !---------------------------------------------------------!
  ! Locals:                                                 !
  ! ys() = array of dependent variables @ stage s of step i !
  ! ki() = array of Runge-Kutta k/h                         !
  ! ci   = Butcher Table c-vector                           !
  ! aij  = Butcher Table A-matrix                           !
  ! bi   = Butcher Table b-vector                           !
  ! dbi  = b-b* vector difference for error calculation     !
  !=========================================================!
  implicit none
  real(dp), intent(in)          :: t, h, cell_volume, velocity_offset(:,:), xi(:)
  integer, intent(in)           :: idx(:)
  type(FluidCell), intent(inout):: grid(:,:,:)
  real(dp), dimension(size(xi)) :: k1, k2, k3, k4, k5, k6, xs
  real(dp)                      :: c2, c3, c4, c5, c6, a21, a31, a32, a41, &
                                   a42, a43, a51, a52, a53, a54, a61, a62, &
                                   a63, a64, a65, b1, b3, b4, b6
  parameter ( c2 = .2_dp, c3 = .3_dp, c4 = .6_dp, c5 = 1._dp, c6 = .875_dp, &
              a21 = .2_dp, a31 = .075_dp, a32 = .225_dp, &
              a41 = .3_dp, a42 = -.9_dp, a43 = 1.2_dp, &
              a51 = -11._dp / 54._dp, a52 = 2.5_dp, a53 = -70._dp / 27._dp, &
              a54 = 35._dp / 27._dp, a61 = 1631._dp / 55296._dp, &
              a62 = 175._dp / 512._dp, a63 = 575._dp / 13824._dp, &
              a64 = 44275. / 110592._dp, a65 = 253. / 4096._dp, &
              b1 = 37._dp / 378._dp, b3 = 250._dp / 621._dp, &
              b4 = 125._dp / 594._dp, b6 = 512._dp / 1771._dp )

  interface interpolated_velocities
    subroutine interp_vel(grid, x, idx, cell_volume, velocity_offset)
      use fluid_types
      implicit none
      type(FluidCell), intent(inout) :: grid(:,:,:) ! Grid
      real(dp), intent(in)           :: x(:), cell_volume, velocity_offset(:,:) ! Coordinates,
      integer, intent(in)            :: idx(:)
    end subroutine interp_vel
  end interface interpolated_velocities

  ! Calculate k1.
  call interp_vel(grid, xi, idx, cell_volume, velocity_offset)
  k1 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  xs = xi + h * a21 * k1

  ! Calculate k2.
  call interp_vel(grid, xs, idx, cell_volume, velocity_offset)
  k2 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  xs = xi + h * (a31 * k1 + a32 * k2)

  ! Calculate k3.
  call interp_vel(grid, xs, idx, cell_volume, velocity_offset)
  k3 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  xs = xi + h * (a41 * k1 + a42 * k2 + a43 * k3)

  ! Calculate k4.
  call interp_vel(grid, xs, idx, cell_volume, velocity_offset)
  k4 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  xs = xi + h * (a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4)

  ! Calculate k5.
  call interp_vel(grid, xs, idx, cell_volume, velocity_offset)
  k5 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  xs = xi + h * (a61 * k1 + a62 * k2 + a63 * k3 + a64 * k4 + a65 * k5)

  ! Calculate k6.
  call interp_vel(grid, xs, idx, cell_volume, velocity_offset)
  k6 = grid(idx(1),idx(2),idx(3)) % tmp_velocity
  grid(idx(1),idx(2),idx(3)) % tmp_x = xi + h * (b1 * k1 + b3 * k3 + b4 * k4 + b6 * k6)

end subroutine rkck



end module fluid
