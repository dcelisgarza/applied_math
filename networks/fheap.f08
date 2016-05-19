module fheap
  private
  public :: FibNode, FibHeap
  ! Fibonacci heap originally written in C++ by beniz.
  ! https://github.com/beniz/fiboheap
  ! OOP fortran linked lists
  ! https://www.pgroup.com/lit/articles/insider/v3n2a2.htm

  type FibNode
    private
    integer                :: key
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

  type, extends(FibNode)    :: FibHeap
    integer                 :: n
    class(FibNode), pointer :: nmin
  end type FibHeap


contains

  recursive subroutine delete_fibnodes(x)
    ! Delete all nodes.
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

  subroutine insert(x, nmin)
    ! Insert new node to fiboheap
    implicit none
    class(FibNode), intent(inout), pointer :: x
    class(FibHeap), intent(inout), pointer :: nmin

    ! Set up the properties of the new node.
    x % degree = 0
    x % parent => null()
    x % child  => null()
    x % mark   = .false.

    ! Check if nMin (node with the minimum key) is Associated.
    ma: if ( .not. associated(nmin) ) then
      ! If nmin is not associated, then there are no nodes in the heap.
      ! Therefore, there is nothing else for x to point to.
      ! We associate x%right and x%left to x (self-reference).
      x % right => x
      x % left  => x % right
      ! We associate nmin with x%left. Through inheritance, it points to x.
      nmin % nmin => x % left
    else ma
      ! Insert x into root list.
      nmin % left % right => x
                x % left  => nmin % left
             nmin % left  => x
              x   % right => nmin % nmin

      if (x % key < nmin % key) nmin % nmin => x
    end if ma

    nmin % n = nmin % n + 1
  end subroutine insert

  subroutine clear_heap(nmin)
    ! Clear all nodes, creating a new heap.
    implicit none
    class(FibHeap), intent(inout), pointer :: nmin

    call delete_fibnodes(nmin % nmin)
    nmin % n = 0
  end subroutine clear_heap

  function get_nmin(nmin)
    ! Get minimum node.
    implicit none
    class(FibHeap), intent(in), pointer :: nmin
    class(FibHeap), pointer             :: get_nmin

    get_nmin => nmin

  end function get_nmin


end module fheap
