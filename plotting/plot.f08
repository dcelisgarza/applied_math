module plot
  implicit none
contains
  subroutine wxtterm(filename, plot_size, font, font_size)
    implicit none
    character(len=*), intent(in):: filename, plot_size, font, font_size

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
end module plot
