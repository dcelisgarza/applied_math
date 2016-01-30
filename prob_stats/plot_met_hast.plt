set terminal pngcairo color enhanced
set output 'plot.png'
set size ratio 1
set xrange [-20:20]
set yrange [-20:20]
plot 'data.dat' w d
