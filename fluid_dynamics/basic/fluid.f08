module fluid
  use nrtype
  use ode
  implicit none

  ! Based on
  ! http://cg.informatik.uni-freiburg.de/intern/seminar/gridFluids_fluid_flow_for_the_rest_of_us.pdf

  ! Define type for a fluid.
  type FluidQty
    ! Value of the quantity.
    real(dp) :: value
    !--------------------------------------------------------------------------!
    ! Offset of the quantity with respect to the main vertex.
    ! offset(1) := offset_x, offset(2) := offset_y, offset(3) := offset_z
    ! Pressure is offset to the middle of the cell:
    ! coordinates( Pressure ) = ( x + d/2,  y + d/2, z + d/2)
    ! Velocity is offset to the center of the minimal surface of its coordinate.
    ! coordinates( Vel x ) = ( x        , y + d/2, z + d/2 )
    ! coordinates( Vel y ) = ( x + d/2, y        , z + d/2 )
    ! coordinates( Vel z ) = ( x + d/2, y + d/2, z         )
    real(dp) :: offset(3)
  end type FluidQty

  type FluidBuffer
    integer :: cell_type
    logical :: initial = .true.
    integer :: layer
    type(FluidQty) :: velocity(3)
    type(FluidQty) :: pressure
  end type FluidBuffer

  type FluidCell
    !--------------------------------------------------------------------------!
    ! When using a hash table, make a type for marker particles.
    ! Marker particle.
    ! marker = 0 -> no particle, marker = 1 -> particle
    integer :: marker
    ! Cell type.
    ! m = 1 -> empty.
    ! m = A_n -> material n, n > 1
    ! In this example, m = 2 -> fluid, m = 3 -> solid
    integer :: cell_type
    !--------------------------------------------------------------------------!
    ! Help simulation step that builds outwards from the fluid into the air buffer
    integer :: layer
    !--------------------------------------------------------------------------!
    ! Coordinates of the main vertex of the cell.
    ! Main vertex := left-down of the 2D cell, or left-down-back of the 3D cell.
    ! x(1) := x, x(2) := y, x(3) := z
    real(dp) :: x(3)
    !--------------------------------------------------------------------------!
    ! Dimensions of the cell.
    real(dp) :: cell_dim
    !--------------------------------------------------------------------------!
    ! Velocity.
    ! velocity(1) = vel_x, velocity(2) = vel_y, velocity(3) = vel_z
    type(FluidQty) :: velocity(3)
    !--------------------------------------------------------------------------!
    ! Pressure.
    type(FluidQty) :: pressure
    !--------------------------------------------------------------------------!
    ! Buffer cell
    type(FluidBuffer) :: buffer(6)
  end type FluidCell

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
          tmp_vel = norm2(grid(i, j, k) % velocity % value)
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

    ! TO DO: MAKE THE PARTICLE MARKERS A LINKED LIST WITH THE INDICES OF THE CELL IN WHICH THEY ARE FOUND. THIS TRIPLE LOOP SHIT CAN GO TO HELL, AND YOU CAN SIMPLY TRAVERSE THE LINKED LIST AND OBTAINING THE INDICES OF THE CORRESPONDING GRID CELL. THIS CAN HELP WHEN EXPANDING THE CODE TO INCLUDE SINKS AND SOURCES.

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

end module fluid
