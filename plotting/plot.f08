module plot
implicit none

  type titles
    character(:), allocatable :: title
    integer                   :: length
  end type titles

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
  character(len=*), intent(in):: filename
  character(len=*), intent(in), optional :: font, plot_name
  integer, intent(in), optional          :: plot_size(2), font_size

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

  if (present(plot_name) .eqv. .true.) write(1,*) "set output '" // plot_name // ".png'"
end subroutine pngterm

subroutine svgterm(filename, plot_name, plot_size, font, font_size)
  implicit none
  character(len=*), intent(in):: filename, plot_name
  character(len=*), intent(in), optional :: font
  integer, intent(in), optional          :: plot_size(2), font_size

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

subroutine plt_labels(xlabel,ylabel,zlabel,title)
  implicit none
  character(len=*), intent(in), optional :: xlabel, ylabel, zlabel, title

  if (present(xlabel) .eqv. .true.) write(1,*) "set xlabel '" // xlabel // "'"
  if (present(ylabel) .eqv. .true.) write(1,*) "set ylabel '" // ylabel // "'"
  if (present(zlabel) .eqv. .true.) write(1,*) "set zlabel '" // zlabel // "'"
  if (present(title)  .eqv. .true.) write(1,*) "set title '"  // title  // "'"
end subroutine plt_labels

subroutine scale(xscale,yscale,zscale)
  implicit none
  character(len=*), intent(in), optional :: xscale, yscale, zscale

  if (present(xscale) .eqv. .true.) write(1,*) "set " // xscale
  if (present(yscale) .eqv. .true.) write(1,*) "set " // yscale
  if (present(zscale) .eqv. .true.) write(1,*) "set " // zscale
end subroutine scale

subroutine range(xrange,yrange,zrange)
  implicit none
  real(4), intent(in), optional :: xrange(2), yrange(2), zrange(2)

  if (present(xrange) .eqv. .true.) write(1,*) "set xrange [", xrange(1), ": ", xrange(2), "]"
  if (present(yrange) .eqv. .true.) write(1,*) "set yrange [", yrange(1), ": ", yrange(2), "]"
  if (present(zrange) .eqv. .true.) write(1,*) "set zrange [", zrange(1), ": ", zrange(2), "]"
end subroutine range

! For axis formats
! fmt
! if (present(fmt)    .eqv. .true.) write(1,*) "set format '" // fmt    // "'"

subroutine ticks(xticks,yticks,zticks,mxticks,myticks,mzticks)
  implicit none
  real(4), intent(in), optional :: xticks(:), yticks(:), zticks(:)
  integer, intent(in), optional :: mxticks, myticks, mzticks

  check_xtick: if (present(xticks) .eqv. .true.) then
    if (size(xticks) > 3 .or. size(xticks) < 1) write(*,*) " Error:: plot.f08:: ticks: check_xtick: &
                                                               1 <= size(xticks) <= 3, current size = ", size(xticks)
    check_xsize: if (size(xticks) == 1) then
      write(1,*) "set xtics " , xticks(1)
    else if (size(xticks) == 2) then check_xsize
      write(1,*) "set xtics " , xticks(1) , ", ", xticks(2)
    else
      write(1,*) "set xtics " , xticks(1) , ", ", xticks(2) , ", ", xticks(3)
    endif check_xsize
  endif check_xtick

  check_ytick: if (present(yticks) .eqv. .true.) then
    if (size(yticks) > 3 .or. size(yticks) < 1) write(*,*) " Error:: plot.f08:: ticks: check_xtick: &
                                                               1 <= size(yticks) <= 3, current size = ", size(yticks)
    check_ysize: if (size(yticks) == 1) then
      write(1,*) "set ytics " , yticks(1)
    else if (size(yticks) == 2) then
      write(1,*) "set ytics " , yticks(1) , ", ", yticks(2)
    else if (size(yticks) == 3) then check_ysize
      write(1,*) "set ytics " , yticks(1) , ", ", yticks(2) , ", ", yticks(3)
    endif check_ysize
  endif check_ytick

  check_ztick: if (present(zticks) .eqv. .true.) then
    if (size(zticks) > 3 .or. size(zticks) < 1) write(*,*) " Error:: plot.f08:: ticks: check_xtick: &
                                                               1 <= size(zticks) <= 3, current size = ", size(zticks)
    check_zsize: if (size(zticks) == 1) then
      write(1,*) "set ztics " , zticks(1)
    else if (size(zticks) == 2) then
      write(1,*) "set ztics " , zticks(1) , ", ", zticks(2)
    else if (size(zticks) == 3) then check_zsize
      write(1,*) "set ztics " , zticks(1) , ", ", zticks(2) , ", ", zticks(3)
    endif check_zsize
  endif check_ztick

  if (present(mxticks) .eqv. .true.) write(1,*) "set mxtics ", mxticks
  if (present(myticks) .eqv. .true.) write(1,*) "set mytics ", myticks
  if (present(mzticks) .eqv. .true.) write(1,*) "set mztics ", mzticks
end subroutine ticks

subroutine grid(xticks,yticks,zticks,linestyle)
  implicit none
  ! Grid ticks := 1D Arrays of size 2.
  ! [major, minor]
  ! 0 = not using tick, 1 = using tick type.
  ! E.g. [0,1] := not using major ticks, using minor ticks; [1,0] := using major ticks, not using minor ticks.
  ! Line style := 1D array of size 6.
  ! [xtick_ls, mxtick_ls, ytick_ls, mytick_ls, ztick_ls, mztick_ls]
  ! 0 = not using line style (automatic)
  integer, intent(in), optional :: xticks(2), yticks(2), zticks(2), linestyle(6)

  check_xticks: if (present(xticks) .eqv. .true.) then
    xtick_ls: if (xticks(1) == 1 .and. linestyle(1) /= 0) then
      write(1,*) "set grid xtics ls ", linestyle(1)
    else if (xticks(1) == 1 .and. linestyle(1) == 0) then xtick_ls
      write(1,*) "set grid xtics"
    end if xtick_ls

    mxtick_ls: if (xticks(2) == 1 .and. linestyle(2) /= 0) then
      write(1,*) "set grid mxtics ls ", linestyle(2)
    else if (xticks(2) == 1 .and. linestyle(2) == 0) then mxtick_ls
      write(1,*) "set grid mxtics"
    end if mxtick_ls
  end if check_xticks

  check_yticks: if (present(yticks) .eqv. .true.) then
    ytick_ls: if (yticks(1) == 1 .and. linestyle(3) /= 0) then
      write(1,*) "set grid ytics ls ", linestyle(3)
    else if (yticks(1) == 1 .and. linestyle(3) == 0) then ytick_ls
      write(1,*) "set grid ytics"
    end if ytick_ls

    mytick_ls: if (yticks(2) == 1 .and. linestyle(4) /= 0) then
      write(1,*) "set grid mytics ls ", linestyle(4)
    else if (yticks(2) == 1 .and. linestyle(4) == 0) then mytick_ls
      write(1,*) "set grid mytics"
    end if mytick_ls
  end if check_yticks

  check_zticks: if (present(zticks) .eqv. .true.) then
    ztick_ls: if (zticks(1) == 1 .and. linestyle(5) /= 0) then
      write(1,*) "set grid ztics ls ", linestyle(5)
    else if (zticks(1) == 1 .and. linestyle(5) == 0) then ztick_ls
      write(1,*) "set grid ztics"
    end if ztick_ls

    mztick_ls: if (zticks(2) == 1 .and. linestyle(6) /= 0) then
      write(1,*) "set grid mztics ls ", linestyle(6)
    else if (zticks(2) == 1 .and. linestyle(6) == 0) then mztick_ls
      write(1,*) "set grid mztics"
    end if mztick_ls
  end if check_zticks
end subroutine grid

subroutine linestyle(l, lc, lt, lw, pt, ps)
  implicit none
  integer, intent(in)                     :: l
  character(len=*), intent(in), optional  :: lc
  integer, intent(in), optional           :: lt, pt
  real(4), intent(in), optional           :: lw, ps

  write(1,*) "set style line", l, " \"
  if (present(lc) .eqv. .true.) write(1,*) "lc " // lc // " \"
  if (present(lt) .eqv. .true.) write(1,*) "lt ", lt, " \"
  if (present(lw) .eqv. .true.) write(1,*) "lw ", lw, " \"
  if (present(pt) .eqv. .true.) write(1,*) "pt ", pt, " \"
  if (present(pt) .eqv. .true.) write(1,*) "ps ", ps
end subroutine linestyle

subroutine dplot2d(filename,using,nplots)
  implicit none
  character(len=*), intent(in)   :: filename
  integer, intent(in), optional  :: using(:), nplots
  integer                        :: i, j
  character(len=len(filename)+4) :: dfilename

  dfilename = filename // ".dat"

  j = 0
  write(1,*) "plot '" // dfilename // "' \"
  columns: if (present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,*) ", '" // dfilename // "' \"
      write(1,*) "u ", using(i+j), ": ", using(i+1+j), " \"
      ! Add optional parameter for lines, linepoints or points.
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
  write(1,*) "splot '" // dfilename // "' \"
  columns: if (present(using) .eqv. .true.) then
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,*) ", '" // dfilename // "' \"
      write(1,*) "u ", using(i+j), ": ", using(i+1+j), ": ", using(i+2+j), " \"
      ! Add an optional linear array of size 2, with one entry being deciding whether to use lines, linepoints or points and the second being the corresponding line style. And don't forget about the title!
      j = j + 2
    end do plots_loop
  end if columns
end subroutine dplot3d

subroutine adplot3d(dfilename,pngname,interval,step,using,nplots,title)
  implicit none
  character(len=*), intent(in)       :: dfilename, pngname
  integer, intent(in)                :: interval(2)
  integer, intent(in), optional      :: step, using(:), nplots
  character(len=*), optional         :: title(:)
  integer                            :: i, j
  character(len=len(dfilename)+4)    :: cdfilename
  ! EXPLANATION OF EVERY http://xmodulo.com/how-to-plot-using-specific-rows-of-data-file-with-gnuplot.html
  ! plot "my.dat" every A:B:C:D:E:F
  ! A: line increment
  ! B: data block increment
  ! C: The first line
  ! D: The first data block
  ! E: The last line
  ! F: The last data block
  cdfilename = dfilename // ".dat"

  write(1,*) "n = 0"
  if (present(step) .eqv. .true.) then
    write(1,*) "do for [i = ", interval(1), ":", interval(2), ":", step, "] {"
  else
    write(1,*) "do for [i = ", interval(1), ":", interval(2), "] {"
  end if
  write(1,*) "n = n + 1"
  write(1,*) "set output sprintf('tmp/"//pngname//"%d.png',n)"

  write(1,*) "splot '" // cdfilename // "' \"
  columns: if (present(using) .eqv. .true.) then
    j = 0
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,*) ", '" // cdfilename // "' \"
      write(1,*) "u ", using(i+j), ": ", using(i+1+j), ": ", using(i+2+j), " \"
      write(1,*) "every ::", interval(1), "::i w l ls ", i," \"
      write(1,*) "notitle \"
      write(1,*) ", '"// cdfilename //"' u ", using(i+j), ": ", using(i+1+j), ": ", using(i+2+j), " \"
      write(1,*) "every ::i::i w p ls ", i," \"
      write(1,*) "title '", trim(title(i)), "' \"
      j = j + 2
    end do plots_loop
  else columns
    write(1,*) "every ::", interval(1), "::i w l ls 1 \"
    write(1,*) ", '"// cdfilename //" every ::i::i w p ls 1 \"
  end if columns
  write(1,*) "}"
end subroutine adplot3d

subroutine adplot2d(dfilename,pngname,interval,step,using,nplots)
  implicit none
  character(len=*), intent(in)   :: dfilename, pngname
  integer, intent(in), optional  :: step, using(:), nplots
  integer, intent(in)            :: interval(2)
  integer                        :: i, j
  character(len=len(dfilename)+4) :: cdfilename
  ! EXPLANATION OF EVERY http://xmodulo.com/how-to-plot-using-specific-rows-of-data-file-with-gnuplot.html
  ! plot "my.dat" every A:B:C:D:E:F
  ! A: line increment
  ! B: data block increment
  ! C: The first line
  ! D: The first data block
  ! E: The last line
  ! F: The last data block
  cdfilename = dfilename // ".dat"

  write(1,*) "n = 0"
  if (present(step) .eqv. .true.) then
    write(1,*) "do for [i = ", interval(1), ":", interval(2), ":", step, "] {"
  else
    write(1,*) "do for [i = ", interval(1), ":", interval(2), "] {"
  end if
  write(1,*) "n = n + 1"
  write(1,*) "set output sprintf('tmp/"//pngname//"%d.png',n)"

  write(1,*) "splot '" // cdfilename // "' \"
  columns: if (present(using) .eqv. .true.) then
    j = 0
    plots_loop: do i = 1, nplots
      if (i > 1) write(1,*) ", '" // cdfilename // "' \"
      write(1,*) "u ", using(i+j), ": ", using(i+1+j), " \"
      write(1,*) "every ::", interval(1), "::i w l ls ", i," \"
      write(1,*) ", '"// cdfilename //"' u ", using(i+j), ": ", using(i+1+j), " \"
      write(1,*) "every ::i::i w p ls ", i," \"
      j = j + 2
    end do plots_loop
  else columns
    write(1,*) "every ::", interval(1), "::i w l ls 1 \"
    write(1,*) ", '"// cdfilename //" every ::i::i w p ls 1 \"
  end if columns
  write(1,*) "}"
end subroutine adplot2d

end module plot
