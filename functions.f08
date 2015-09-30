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
end module functions
