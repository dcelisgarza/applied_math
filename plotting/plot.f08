module plot
implicit none
contains

! GNUPlot Terminals.
subroutine int1_char(cint,int)
  implicit none
  integer, intent(in) :: int
  character(1)        :: cint

  ! This is how you calculate the number of integers and decimals in the format statement.
  !int = nint( maxval(plot_size) )
  !dec = nint( mod( maxval(plot_size), 1. ) / 0.001 )

  int_len: if (1 > int) then
    cint = "0"
  else if (1 <= int .and. int < 10) then int_len
    cint = "1"
  else if (10 <= int .and. int < 100) then int_len
    cint = "2"
  else if (100 <= int .and. int < 1000) then int_len
    cint = "3"
  end if int_len
end subroutine int1_char

subroutine wxtterm(filename, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in)           :: filename
  character(len=*), intent(in), optional :: font
  integer, intent(in), optional          :: plot_size(2)
  integer, intent(in), optional          :: font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal wxt \"

  plt_size: if (present(plot_size) .eqv. .true.) then
    write(1,*) "size ", plot_size(1), ", " , plot_size(2), " \"
  end if plt_size

  write(1,*) "enhanced font \"

  fnt: if (present(font) .eqv. .true.) then
    write(1,*) "'" // font, ' \'
    fnt_size: if (present(font_size) .eqv. .true.) then
      write(1,*) ", ", font_size, " \"
    end if fnt_size
    write(1,*) "' \"
  end if fnt

   write(1,*) "persist"

end subroutine wxtterm

subroutine pngterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name
  character(len=*), intent(in), optional :: font
  integer, intent(in), optional          :: plot_size(2)
  integer, intent(in), optional          :: font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal pngcairo \"

  plt_size: if (present(plot_size) .eqv. .true.) then
    write(1,*) "size ", plot_size(1), ", " , plot_size(2), " \"
  end if plt_size

  write(1,*) "enhanced font \"

  fnt: if (present(font) .eqv. .true.) then
    write(1,*) "'" // font, ' \'
    fnt_size: if (present(font_size) .eqv. .true.) then
      write(1,*) ", ", font_size, " \"
    end if fnt_size
    write(1,*) "'"
  end if fnt

  write(1,*) "set output '" // plot_name // ".png'"
end subroutine pngterm

subroutine svgterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name
  character(len=*), intent(in), optional :: font
  integer, intent(in), optional          :: plot_size(2)
  integer, intent(in), optional          :: font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal svg \"

  plt_size: if (present(plot_size) .eqv. .true.) then
    write(1,*) "size ", plot_size(1), ", " , plot_size(2), " \"
  end if plt_size

  write(1,*) "enhanced font \"

  fnt: if (present(font) .eqv. .true.) then
    write(1,*) "'" // font, ' \'
    fnt_size: if (present(font_size) .eqv. .true.) then
      write(1,*) ", ", font_size, " \"
    end if fnt_size
    write(1,*) "'"
  end if fnt

  write(1,*) "set output '" // plot_name // ".svg'"
end subroutine svgterm

subroutine epsterm(filename, plot_name, plot_size, plot_size_units, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name
  character(len=*), intent(in), optional :: font, plot_size_units
  real(4), intent(in), optional          :: plot_size(2)
  integer, intent(in), optional          :: font_size

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal postscript eps \"

  plt_size: if (present(plot_size) .eqv. .true.) then
    write(1,*) "size ", plot_size(1), " \"
    plt_size_units: if (present(plot_size_units) .eqv. .true.) then
      write(1,*) plot_size_units // ", ", plot_size(2), plot_size_units, " \"
    else plt_size_units
      write(1,*) ", ", plot_size(2), " \"
    end if plt_size_units
  end if plt_size

  write(1,*) "enhanced color \"

  fnt: if (present(font) .eqv. .true.) then
    write(1,*) "font '" // font, ' \'
    fnt_size: if (present(font_size) .eqv. .true.) then
      write(1,*) ", ", font_size, " \"
    end if fnt_size
    write(1,*) "' \"
  end if fnt

  write(1,*) "set output '" // plot_name // ".eps'"
end subroutine epsterm

subroutine epslatexterm(filename, plot_name, plot_size, plot_size_units)
  implicit none
  character(len=*), intent(in):: filename, plot_name
  character(len=*), intent(in), optional :: plot_size_units
  real(4), intent(in), optional          :: plot_size(2)

  open(unit = 1, file = filename//'.gnu')
  write(1,*) "set terminal epslatex \"

  plt_size: if (present(plot_size) .eqv. .true.) then
    write(1,*) "size ", plot_size(1), " \"
    plt_size_units: if (present(plot_size_units) .eqv. .true.) then
      write(1,*) plot_size_units // ", ", plot_size(2), plot_size_units, " \"
    else plt_size_units
      write(1,*) ", ", plot_size(2), " \"
    end if plt_size_units
  end if plt_size

  write(1,*) "color colortext"

  write(1,*) "set output '" // plot_name // ".tex'"
end subroutine epslatexterm

subroutine format(xlabel,ylabel,zlabel,title,fmt)
  implicit none
  character(len=*), intent(in), optional :: xlabel, ylabel, zlabel, title, fmt

  if (present(xlabel) .eqv. .true.) write(1,*) "set xlabel '" // xlabel // "'"
  if (present(ylabel) .eqv. .true.) write(1,*) "set ylabel '" // ylabel // "'"
  if (present(zlabel) .eqv. .true.) write(1,*) "set xlabel '" // zlabel // "'"
  if (present(title)  .eqv. .true.) write(1,*) "set title '"  // title  // "'"
  if (present(fmt)    .eqv. .true.) write(1,*) "set format '" // fmt    // "'"
end subroutine format

subroutine ticks(xticks,yticks,zticks)
  implicit none
  real(4), intent(in), optional :: xticks(3), yticks(3), zticks(3)

  if (present(xticks) .eqv. .true.) write(1,*) "set xtics ", xticks(1), ", ", xticks(2), ", ", xticks(3)
  if (present(yticks) .eqv. .true.) write(1,*) "set ytics ", yticks(1), ", ", yticks(2), ", ", yticks(3)
  if (present(zticks) .eqv. .true.) write(1,*) "set xtics ", zticks(1), ", ", zticks(2), ", ", zticks(3)
end subroutine ticks

subroutine range(xrange,yrange,zrange)
  implicit none
  real(4), intent(in), optional :: xrange(2), yrange(2), zrange(2)

  if (present(xrange) .eqv. .true.) write(1,*) "set xrange [", xrange(1), ": ", xrange(2), "]"
  if (present(yrange) .eqv. .true.) write(1,*) "set yrange [", yrange(1), ": ", xrange(2), "]"
  if (present(zrange) .eqv. .true.) write(1,*) "set xrange [", zrange(1), ": ", xrange(2), "]"
end subroutine range

subroutine scale(xscale,yscale,zscale)
  implicit none
  character(len=*), intent(in), optional :: xscale, yscale, zscale

  if (present(xscale) .eqv. .true.) write(1,*) "set " // xscale
  if (present(yscale) .eqv. .true.) write(1,*) "set " // yscale
  if (present(zscale) .eqv. .true.) write(1,*) "set " // zscale
end subroutine scale

subroutine linestyle(l, lc, lt, lw, pt, ps)
  implicit none
  integer, intent(in)                     :: l
  character(len=*), intent(in), optional  :: lc
  integer, intent(in), optional           :: lt, pt
  real(4), intent(in), optional           :: lw, ps

  write(1, "(A,i2)", advance = "no") " set style line", l, " "
  if (present(lc) .eqv. .true.) write(1, "(A)",     advance = "no") "lc " // lc
  if (present(lt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "lt ", lt
  if (present(lw) .eqv. .true.) write(1, "(A, i2)", advance = "no") "lw ", lw
  if (present(pt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "pt ", pt
  if (present(pt) .eqv. .true.) write(1, "(A, i2)", advance = "no") "ps ", ps
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
  columns: if (present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,"(A)", advance = "no") ", '" // dfilename // "' "
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
  columns: if (present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,"(A)", advance = "no") ", '" // dfilename // "' "
      write(1,"(3(A, i3))", advance = "no") "u ", using(i+j), ": ", using(i+1+j), ": ", using(i+2+j)
      j = j + 2
    end do plots_loop
  end if columns
end subroutine dplot3d

subroutine danim2d(dfilename,interval,pngname)
  implicit none
  character(len=*), intent(in)   :: dfilename, pngname ! Data file name, magnitu
  integer, intent(in), optional  :: interval(2)
  integer                        :: i, j, k
  character(len=len(dfilename)+4) :: cdfilename ! Complete data filename

  ! Getting the data file name
  cdfilename = dfilename // ".dat"

  write(1,*) "n = 0"
  write(1,*) "do for [i = ", interval(1)-1, ":", interval(2), "] {"
  write(1,*) "n = n + 1"
  write(1,*) "set output sprintf('tmp/"//pngname//"%d.png',n)"
  ! EXPLANATION OF EVERY http://xmodulo.com/how-to-plot-using-specific-rows-of-data-file-with-gnuplot.html
  ! plot "my.dat" every A:B:C:D:E:F
  ! A: line increment
  ! B: data block increment
  ! C: The first line
  ! D: The first data block
  ! E: The last line
  ! F: The last data block
  write(1,*) "plot '"//cdfilename//"' every ::", interval(1)-1, "::i w l, \"
  ! Condition this later on to make it visible under user command.
  write(1,*) "'"//cdfilename//"' every ::i::i w p"
  write(1,*) "}"

end subroutine danim2d

subroutine danim3d(filename,interval,mag_frames,pngname)
  implicit none
  character(len=*), intent(in)   :: filename, mag_frames, pngname
  integer, intent(in), optional  :: interval(2)
  integer                        :: i, j, k
  character(len=len(filename)+4) :: dfilename

  ! Getting the data file name
  dfilename = filename // ".dat"

  write(1,*) "n = 0"
  write(1,*) "do for [i = ", interval(1), ":", interval(2), "] {"
  write(1,*) "n = n + 1"
  write(1,*) "set output sprintf('tmp/"//pngname//"%d.png',n)"
  write(1,*) "splot '"//dfilename//"' every :3:", interval(1), "::i w l \"
  ! Condition this later on to make it visible under user command.
  write(1,*) "splot '"//dfilename//"' every :3:i::i w p"
  write(1,*) "}"
end subroutine danim3d

end module plot
