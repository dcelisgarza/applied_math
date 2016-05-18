module adtypes

  type fheap
    ! Fibonacci heap double precision.
    class(*), pointer    :: value => null()
    type(fheap), pointer :: left  => null()
    type(fheap), pointer :: right => null()
    type(fheap), pointer :: up    => null()
    type(fheap), pointer :: down  => null()

  contains
    procedure :: adj  => adj  ! Access adjacent links
    procedure :: fmin => fmin ! Find minimum.
    procedure :: join => join ! Join.
    procedure :: ins  => ins  ! Insert.
    procedure :: dmin => dmin ! Delete.
    procedure :: dkey => dkey ! Decrease key.
    procedure :: del  => del  ! Delete node.
  end type fheap

contains

  function adj(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: adj
    adj % left  => this % left
    adj % right => this % right
    adj % up    => this % up
    adj % down  => this % down
  end function adj

  function fmin(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: fmin
  end function fmin

  function join(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: join
    !join% => this%
  end function join

  function ins(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: ins
  end function ins

  function dmin(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: dmin
  end function dmin

  function dkey(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: dkey
  end function dkey

  function del(this)
    implicit none
    class(fheap) :: this
    class(fheap), pointer :: del
  end function del

end module adtypes
