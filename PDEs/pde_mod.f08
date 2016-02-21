module finite_elements
  use nrtype
contains

  subroutine fepdelin(u,ncyc)
    implicit none
    real(dp), intent(inout) :: u(:,:) ! Function u.
    integer, intent(in)     :: ncyc   ! Number of cycles.
    integer                 :: i, n, ng, cn, cng ! Position, grid dimension.
    type a2d
      real(dp), allocatable :: a(:,:)
    end type a2d
    type(a2d), allocatable  :: rho(:)
    real(dp), dimension(:,:), allocatable :: ui, uj

    ! Exit if both dimensions aren't equal. We need a square matrix.
    dimu: if(size(u,1) /= size(u,2)) then
      write(*,*) ' Error: subr:: fepdelin:: dimu. U must be an N x N matrix.'
      return
    end if dimu

    ! Number of multigrids.
    ng = nint(log(n-1.)/log(2.))

    ! Exit if n-1 isn't a power of two.
    dimn: if (n /= 2**ng + 1) then
      write(*,*) ' Error: subr:: fepdelin:: dimn. N-1 must be a power of 2.'
    end if dimn

    ! Allocate the rho array.
    allocate( rho(ng) )
    allocate( rho(ng) % a(n,n) )

    ! Assign the array the values of the initial u at the finest grid.
    rho(ng) % a = u

    ! Assigning coarser n and coarser ng.
    cn  = n
    cng = ng
    grid: do
      if (cn <= 3) exit grid
      cn  = cn/2 + 1
      cng = cng - 1
      allocate( rho(cng) % a(cn,cn) )
      rho(cng) % a = hwrsct( rho(cng+1) % p )
    end do grid
  end subroutine fepdelin

  function hwrsct(u)
    ! Function for restricting u to a coarser grid.
    implicit none
    real(dp), intent(in) :: u(:,:)
    real(dp)             :: hwrsct((size(u,1)+1)/2, (size(u,2)+1)/2) ! Restricted u
    integer              :: n, cn ! n and coarser n

    ! Current grid dimension.
    n = size(u,1)
    ! Coarser grid dimension.
    cn = (n+1)/2

    ! Boundary conditions.
    hwrsct(1:cn,1)  = u(1:n:2,1)
    hwrsct(1:cn,cn) = u(1:n:2,n)
    hwrsct(1,1:cn)  = u(1,1:n:2)
    hwrsct(cn,1:cn) = u(n,1:n:2)

    ! Internal elements.
    hwrsct(2:cn-1,2:cn-1) = 0.5_dp * u(3:n-2:2, 3:n-2:2) + 0.125_dp * &
                            (&
                            u(4:n-1:2,3:n-1:2)&
                            )
  end function hwrsct

end module finite_elements
