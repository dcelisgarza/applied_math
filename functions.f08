module functions
implicit none
integer, parameter :: dp = kind(1.d0)
contains
  function heaviside(x)
    !=======================================================!
    ! Heaviside function                                    !
    ! Daniel Celis Garza 15 Jan 2015                        !
    !=======================================================!
    real(dp) :: heaviside, x
    heaviside = 0.0_dp
    if (x .ge. 0.0) heaviside = 1.0_dp
  end function heaviside

  function kdelta(i,j)
    !=======================================================!
    ! Kroencker delta function                              !
    ! Daniel Celis Garza 15 Jan 2015                        !
    ! delta_ij                                              !
    !=======================================================!
    integer i, j
    real(dp) kdelta
    kdelta = 0.0_dp
    if (i .eq. j) kdelta = 1.0_dp
  end function kdelta
end module functions
