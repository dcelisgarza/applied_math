module plot
implicit none
contains

! GNUPlot Terminals.
subroutine wxtterm(filename, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in) :: filename, plot_size, font, font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal wxt size " // plot_size // " enhanced font '" // font // "," // font_size // "' persist"
end subroutine wxtterm

subroutine pngterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name, plot_size, font, font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal pngcairo size " // plot_size // " enhanced font '" // font // "," // font_size // "'"
  write(1,*) "set output '" // plot_name // ".png'"
end subroutine pngterm

subroutine svgterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name, plot_size, font, font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal svg size " // plot_size // " fname '" // font // " fsize " // font_size
  write(1,*) "set output '" // plot_name // ".svg'"
end subroutine svgterm

subroutine epsterm(filename, plot_name, plot_size, font, font_size, linewidth)
  implicit none
  character(len=*), intent(in):: filename, plot_name, plot_size, font, font_size, linewidth

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal postscript eps size " // plot_size // " enhanced color font '" // font // "," // font_size // "'" &
  // " linewidth " // linewidth
  write(1,*) "set output '" // plot_name // ".eps'"
end subroutine epsterm

subroutine epslatexterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name, plot_size, font, font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal epslatex size " // plot_size // " color colortext " // font_size
  write(1,*) "set output '" // plot_name // ".tex'"
end subroutine epslatexterm

subroutine format(xlabel,ylabel,zlabel,title,fmt)
  implicit none
  character(len=*), intent(in), optional :: xlabel, ylabel, zlabel, title, fmt

  if(present(xlabel) .eqv. .true.) write(1,*) "set xlabel '" // xlabel // "'"
  if(present(ylabel) .eqv. .true.) write(1,*) "set ylabel '" // ylabel // "'"
  if(present(zlabel) .eqv. .true.) write(1,*) "set xlabel '" // zlabel // "'"
  if(present(title)  .eqv. .true.) write(1,*) "set title '"  // title  // "'"
  if(present(fmt)    .eqv. .true.) write(1,*) "set format '" // fmt    // "'"
end subroutine format

subroutine ticks(xticks,yticks,zticks)
  implicit none
  real(4), intent(in), optional :: xticks(3), yticks(3), zticks(3)

  if(present(xticks) .eqv. .true.) write(1,*) "set xtics ", xticks(1), ", ", xticks(2), ", ", xticks(3)
  if(present(yticks) .eqv. .true.) write(1,*) "set ytics ", yticks(1), ", ", yticks(2), ", ", yticks(3)
  if(present(zticks) .eqv. .true.) write(1,*) "set xtics ", zticks(1), ", ", zticks(2), ", ", zticks(3)
end subroutine ticks

subroutine range(xrange,yrange,zrange)
  implicit none
  real(4), intent(in), optional :: xrange(2), yrange(2), zrange(2)

  if(present(xrange) .eqv. .true.) write(1,*) "set xrange [", xrange(1), ": ", xrange(2), "]"
  if(present(yrange) .eqv. .true.) write(1,*) "set yrange [", yrange(1), ": ", xrange(2), "]"
  if(present(zrange) .eqv. .true.) write(1,*) "set xrange [", zrange(1), ": ", xrange(2), "]"
end subroutine range

subroutine scale(xscale,yscale,zscale)
  implicit none
  character(len=*), intent(in), optional :: xscale, yscale, zscale

  if(present(xscale) .eqv. .true.) write(1,*) "set " // xscale
  if(present(yscale) .eqv. .true.) write(1,*) "set " // yscale
  if(present(zscale) .eqv. .true.) write(1,*) "set " // zscale
end subroutine scale

subroutine linestyle(l, lc, lt, lw, pt, ps)
  implicit none
  integer, intent(in)                     :: l
  character(len=*), intent(in), optional  :: lc
  integer, intent(in), optional           :: lt, pt
  real(4), intent(in), optional           :: lw, ps

  write(1, "(A,i2)", advance = "no") " set style line", l, " "
  if(present(lc) .eqv. .true.) write(1, "(A)",     advance = "no") "lc " // lc
  if(present(lt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "lt ", lt
  if(present(lw) .eqv. .true.) write(1, "(A, i2)", advance = "no") "lw ", lw
  if(present(pt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "pt ", pt
  if(present(pt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "ps ", ps
end subroutine linestyle

subroutine dplot2d(filename,using,nplots)
  implicit none
  character(len=*), intent(in)   :: filename
  integer, intent(in), optional  :: using(:), nplots
  integer                        :: i, j
  character(len=len(filename)+4) :: dfilename

  dfilename = filename // ".dat"

  j = 0
  write(1,"(A)", advance = "no") " plot '" // dfilename // "' "
  columns: if(present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if(i > 1) write(1,"(A)", advance = "no") ", '" // dfilename // "' "
      write(1,"(2(A, i3))", advance = "no") "u ", using(i+j), ": ", using(i+1+j)
      j = j + 1
    end do plots_loop
  end if columns
end subroutine dplot2d

subroutine dplot3d(filename,using,nplots)
  implicit none
  character(len=*), intent(in)   :: filename
  integer, intent(in), optional  :: using(:), nplots
  integer                        :: i, j
  character(len=len(filename)+4) :: dfilename

  dfilename = filename // ".dat"

  j = 0
  write(1,"(A)", advance = "no") " splot '" // dfilename // "' "
  columns: if(present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if(i > 1) write(1,"(A)", advance = "no") ", '" // dfilename // "' "
      write(1,"(3(A, i3))", advance = "no") "u ", using(i+j), ": ", using(i+1+j), ": ", using(i+2+j)
      j = j + 2
    end do plots_loop
  end if columns
end subroutine dplot3d

end module plot
