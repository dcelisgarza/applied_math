program plot_main
  use plot
  implicit none
  character(:), allocatable :: filename, plot_name, plot_size, font, font_size, linewidth
  filename = 'test'
  plot_size = '10,2'
  font = 'Helvetica'
  font_size = '14'
  plot_name = 'test'
  linewidth = '2'
  call wxtterm(filename, plot_size, font, font_size)
  call pngterm(filename, plot_name, plot_size, font, font_size)
  call svgterm(filename, plot_name, plot_size, font, font_size)
  call epsterm(filename, plot_name, plot_size, font, font_size, linewidth)
  call epslatexterm(filename, plot_name, plot_size, font, font_size)
end program plot_main
