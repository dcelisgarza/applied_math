module fluid_types
  use nrtype
  ! Define type for a fluid.
  type GridParameters
    ! Offset of the quantity with respect to the main vertex of a given grid cell.
    ! offset(1) := offset_x, offset(2) := offset_y, offset(3) := offset_z
    ! Pressure is offset to the middle of the cell:
    ! coordinates( Pressure ) = ( x + d/2,  y + d/2, z + d/2)
    ! Velocity is offset to the center of the minimal surface of its coordinate.
    ! velocity_offset(1,1:3) := offset( Vel x ) = ( x      , y + d/2, z + d/2 )
    ! velocity_offset(2,1:3) := offset( Vel y ) = ( x + d/2, y      , z + d/2 )
    ! velocity_offset(3,1:3) := offset( Vel z ) = ( x + d/2, y + d/2, z       )
    real(dp) :: velocity_offset(3,3)
    real(dp) :: pressure_offset(3)
    real(dp) :: cell_volume
  end type GridParameters

  type FluidBuffer
    integer :: cell_type
    logical :: initial = .true.
    integer :: layer
    real(dp):: velocity(3)
    real(dp):: pressure
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
    ! Temporary value of position.
    real(dp) :: tmp_x(3)
    !--------------------------------------------------------------------------!
    ! Dimensions of the cell.
    real(dp) :: cell_dim
    !--------------------------------------------------------------------------!
    ! Velocity.
    ! velocity(1) = vel_x, velocity(2) = vel_y, velocity(3) = vel_z
    real(dp) :: velocity(3)
    ! TMP Velocity, advancing velocity fields can't be done in-place.
    real(dp) :: tmp_velocity(3)
    !--------------------------------------------------------------------------!
    ! Pressure.
    real(dp) :: pressure
    !--------------------------------------------------------------------------!
    ! Buffer cell
    type(FluidBuffer) :: buffer(6)
  end type FluidCell
end module fluid_types
