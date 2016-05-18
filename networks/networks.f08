module networks
  use nrtype
  use ieee_arithmetic
!https://en.wikipedia.org/wiki/Dijkstra's_algorithm
!https://en.wikipedia.org/wiki/Brodal_queue
!https://en.wikipedia.org/wiki/Priority_queue
!https://www.cs.auckland.ac.nz/software/AlgAnim/dijkstra.html
!http://www.geeksforgeeks.org/greedy-algorithms-set-6-dijkstras-shortest-path-algorithm/
!http://www.personal.kent.edu/~rmuhamma/Algorithms/MyAlgorithms/GraphAlgor/dijkstraAlgor.htm
!https://www.cs.princeton.edu/~rs/AlgsDS07/15ShortestPaths.pdf
!http://math.mit.edu/~rothvoss/18.304.3PM/Presentations/1-Melissa.pdf
!http://rosettacode.org/wiki/Dijkstra's_algorithm

contains

  function inf(sign)
    implicit none
    integer(i1), intent(in) :: sign
    real                    :: inf

    ! Choose Sign.
    cs: if (sign > 0) then
      inf = ieee_value(inf,  ieee_positive_inf)
    else cs
      inf = ieee_value(inf,  ieee_negative_inf)
    end if cs
  end function inf

  subroutine dijkstra(graph, s, dist, prev)
    implicit none

    integer, intent(in)  :: s ! s = source
    real(dp), intent(in) :: graph(:,:)
    real(dp), intent(out), dimension(size(graph, dim = 1)) :: dist
    integer, dimension(size(graph, dim = 1))               :: visited, prev
    integer :: nv, v, u, ndist, mdist ! v = target node, u = current node with minimum distance to v, ndist = new distance from u to v.
    integer :: i ! counter

    ! Number of vertices (nodes)
    nv = size(graph, dim = 1)

    ! Check Graph Size.
    if ( nv /= size(graph, dim=2) ) &
      write(*,*) " Matrix size error: networks.f08: dijkstra: cgs: Graph provided must be a square matrix. &
      Dimension 1 = ", size(graph, dim=1), " Dimension 2 = ", size(graph, dim=2)

    ! Initialising variables.
    dist = graph(s,1:nv) ! One step distance from node s to the rest.
    !print*, "---------------"
    !print*, dist
    !print*, "---------------"
    prev = -1 ! Previously visited node.
    visited = 0 ! No nodes have been visited.
    visited(s) = 1 ! Source has been visited.
    !print*, "---------------"
    !print*, visited
    !print*, "---------------"

    ! Dijkstra's Algorithm
    da: do i = 1, nv
      if (i == s) cycle da
      call find_nearest(nv, dist, visited, v)
      if (v == -1) return
      visited(v) = 1
      call update_dist(nv, graph, dist, visited, v, prev)

    end do da
  end subroutine dijkstra

  subroutine find_nearest(nv, dist, visited, v)
    implicit none
    integer, intent(in)  :: nv
    real(dp), intent(in) :: dist(:)
    integer, intent(in)  :: visited(:)
    integer, intent(out) :: v
    real(dp) :: mdist
    integer  :: i ! counter

    mdist = inf(1_i1)
    v     = -1

    ! Find Minimum Distance.

    fmd: do i = 1, nv
      ! Check if the node has Not been Visited, and whether it is Connected to the Source node.
      nvc: if (visited(i) == 0 .and. dist(i) < mdist) then
        ! If the node has not been visited and it is connected to the source.
        ! The minimum distance is the distance to node i, and the node is node i.
        mdist = dist(i)
        v     = i
      end if nvc
    end do fmd

  end subroutine find_nearest

  subroutine update_dist(nv, graph, dist, visited, v, prev)
    implicit none
    integer, intent(in)     :: nv
    real(dp), intent(in)    :: graph(:,:)
    integer, intent(in)     :: visited(:)
    integer, intent(in)     :: v
    real(dp), intent(inout) :: dist(:)
    integer, intent(out) :: prev(:)
    real(dp) :: newdist
    integer :: i

    ! Calculate the Minimum Distance.
    cmd: do i = 1, nv
      ! Check if we've visited node i or if node v and i are not connected, move on to the next node i.
      if ( visited(i) == 1 .or. graph(v,i) == inf(1_i1) ) cycle cmd
      ! Check whether the New Distance from node v node i is Less than the Previous minimum distance to node i.
      newdist = dist(v) + graph(v,i)
      ndlpd: if (newdist < dist(i)) then
        ! If it is, then update the new minimum distance.
        dist(i) = newdist
        prev(i) = v
      end if ndlpd
    end do cmd
  end subroutine update_dist

  subroutine print_target(file, t, dist, prev)
    implicit none
    integer, intent(in)  :: file, t, prev(:) ! file unit, target, previous node
    real(dp), intent(in) :: dist(:)

    write(1,*) " Minimum distance = ", dist(t)

    do
      write(1,*) t
      t = prev(t)
      if (t == -1) return
    end do

  end subroutine print_target
end module networks
