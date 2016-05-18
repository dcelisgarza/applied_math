module fheap
  private
  public :: FibNode
  ! Fibonacci heap originally written in C++ by beniz.
  ! https://github.com/beniz/fiboheap
  ! OOP fortran linked lists
  ! https://www.pgroup.com/lit/articles/insider/v3n2a2.htm

  type FibNode
    private
    integer                :: degree
    logical                :: mark
    class(*),      pointer :: value  => null()
    type(FibNode), pointer :: parent => null()
    type(FibNode), pointer :: child  => null()
    type(FibNode), pointer :: left   => null()
    type(FibNode), pointer :: right  => null()

  contains
    procedure :: delete_fibnodes => delete_fibnodes ! Delete all nodes.
  end type FibNode

contains

  recursive subroutine delete_fibnodes(x)
    implicit none
    class(FibNode), intent(inout), target :: x
    class(FibNode), pointer :: cur, tmp

    ! If x is not associated, return. No need to do anything.
    if(.not. associated(x % value)) return

    ! Associate working variable to x.
    cur => x

    ! Loop to Dissociate All Nodes.
    dan: do
      ! Check if Cur % Left is Associated with X.
      clax: if ( .not. associated(cur % left, target = x) ) then
        ! Associate tmp to cur.
        tmp => cur
        ! Associate current to cur % left.
        cur => cur % left
        ! Check whether tmp % child is associated, if it is,
        ! recursively call the function to dissociate it.
        if (associated(tmp % child)) call delete_fibnodes(tmp % child)
        ! Nullify tmp.
        tmp => null()
      else clax
        if (associated(cur % child)) call delete_fibnodes(cur % child)
        cur => null()
        return
      end if clax
    end do dan
  end subroutine delete_fibnodes

end module fheap
