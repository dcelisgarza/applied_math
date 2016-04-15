module lin_alg
  use nrtype
  ! Module for numerical linear algebra.
contains
  subroutine tridiag(c, a, b, y, x)
    implicit none
    ! Tridiagonal matrix
    !
    ! a_(1,1)   a_(2,1)   0         0         ...         0           |
    ! a_(2,1)   a_(2,2)   a_(2,3)   0         ...         0           |
    ! 0         a_(3,2)   a_(3,3)   a_(3,4)   ...         0           |
    ! 0         0         a_(4,3)   a_(4,4)   ...         0
    ! ...       ...       ...       ...       ...         a_(n-1,n)
    ! 0         0         0         ...       a_(n,n-1)   a_(n,n)
    !
    ! Inputs:
    ! c(:) = Lower Diagonal.
    ! a(:) = Main Diagonal.
    ! b(:) = Upper Diagonal.
    ! y(:) = Image.
    ! Outputs:
    ! x(:) = Solution vector.
    ! Locals:
    ! i, j = Counters.
    ! w(:) = Dummy Vector.
    ! c    = Dummy Variable.
    real(dp), intent(in)  :: c(:), a(:), b(:), y(:)
    real(dp), intent(out) :: x(:)
    real(dp) :: w(size(a)), den
    integer  :: n, i

    ! Set the value of n to be the number of entries in the main diagonal.
    n = size(a)
    ! Check that the sizes of the lower, main and upper diagonals make sense.
    check_size: if( n /= size(c) + 1 .or. &
                    n /= size(b) + 1 .or. &
                    n /= size(y)            ) then
      write(*,*) " Error: lin_alg: tridiag: check_size: Sizes of the upper and lower diagonals must be one smaller than the main &
                   diagonal; size(lower_diagonal) = ", size(c), "; size(upper_diagonal) = ", size(b), &
                   "; size(main_diagonal) = ", size(a)
    end if check_size

    ! Calculate the First Denominator.
    den = a(1)
    fd: if ( den == 0._dp ) then
      write(*,*) " Error: lin_alg: tridiag: Next iteration will divide by zero. den = ", den
      stop
    end if fd
    ! Calculate preliminary value of x(1).
    x(1) = y(1) / den

    ! Decomposition and Forward Substitution.
    dfs: do i = 2, n
      w(i) = b(i - 1) / den
      den  = a(i) - c(i - 1) * w(i)
      ! Error if Next Iteration will Divide by 0.
      nid0: if ( den == 0._dp ) then
        write(*,*) " Error: lin_alg: tridiag: dfs: Next iteration will divide by zero. den = ", den
        stop
      end if nid0
      ! x(2) = ( y(2) - c(1) * x(1) ) /
      ! Solution vector.
      x(i) = ( y(i) - c(i - 1) * x(i - 1) ) / den
    end do dfs

    ! Back Substitution.
    bs: do i = n - 1, 1, -1
      x(i) = x(i) - w(i + 1) * x(i + 1)
    end do bs

  end subroutine tridiag
end module lin_alg
