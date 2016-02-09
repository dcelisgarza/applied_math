module nrtype
  integer, parameter :: i1  = selected_int_kind(2)  ! 1-byte integer.
  integer, parameter :: i2  = selected_int_kind(4)  ! 2-byte integer.
  integer, parameter :: i4  = selected_int_kind(9)  ! 4-byte integer.
  integer, parameter :: i8  = selected_int_kind(18) ! 8-byte integer.
  integer, parameter :: i16 = selected_int_kind(22) ! 16-byte integer.
  integer, parameter :: sp  = kind(1.0)           ! Single precision real double.
  integer, parameter :: spc = kind((1.0,1.0))     ! Single precision complex double.
  integer, parameter :: dp  = kind(1.0d0)         ! Max precision real double.
  integer, parameter :: dpc = kind((1.0d0,1.0d0)) ! Max precision complex double.
  integer, parameter :: lgt = kind(.true.)        ! Logical true.
  integer, parameter :: lgf = kind(.false.)       ! Logical false.

  real(dp), parameter :: pi    = 3.141592653589793238462643383279502884197_dp
  real(dp), parameter :: pi2   = 1.57079632679489661923132169163975144209858_dp
  real(dp), parameter :: tau   = 6.283185307179586476925286766559005768394_dp
  real(dp), parameter :: sqrt2 = 1.41421356237309504880168872420969807856967_dp
  real(dp), parameter :: euler = 2.71828182845904523536028747135266249775725_dp
end module nrtype
