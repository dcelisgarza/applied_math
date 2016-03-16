program plot_main
  use plot
  implicit none
  character(:), allocatable :: filename, plot_name, plot_size, font, font_size, linewidth, &
                               xlabel,ylabel,zlabel,title,fmt,xticks,yticks,zticks,xrange,yrange,zrange
  integer :: nplot_size(2)
  integer :: using(6), nplots
  filename = 'test'
  plot_size = '10,2'
  font = 'Helvetica'
  font_size = '14'
  plot_name = 'test'
  linewidth = '2'
  xlabel = 'X-axis'
  ylabel = 'Y-axis'
  zlabel = 'Z-axis'
  title = 'title'
  fmt = 'format'
  xticks = '0,1,0.1'
  yticks = '0,1,0.1'
  zticks = '0,1,0.1'
  xrange = '0:1'
  yrange = '0:1'
  zrange = '0:1'
  using  = [1,2,3,4,5,6]
  nplots = 3
  nplot_size = [500,500]

  call wxtterm(filename, nplot_size, font, 12)
  call pngterm(filename, plot_name, nplot_size, font, 12)
  call svgterm(filename, plot_name, nplot_size, font, 12)
  call epsterm(filename, plot_name, [5.,5.], "cm",font, 12, 3.)
  !call epslatexterm(filename, plot_name, plot_size, font, font_size)
  !call dplot2d(filename,using,nplots)
  !call format(xlabel,ylabel,zlabel,title,fmt)
  !call ticks(xticks,yticks,zticks)
  !call range(xrange,yrange,zrange)

end program plot_main
