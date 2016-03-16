 set terminal wxt \
 size          500 ,          500  \
 enhanced font \
 'Helvetica \
 ,           12  \
 ' \
 persist
 set terminal pngcairo \
 size          500 ,          500  \
 enhanced font \
 'Helvetica \
 ,           12  \
 '
 set output 'test.png'
 set terminal svg \
 size          500 ,          500  \
 enhanced font \
 'Helvetica \
 ,           12  \
 '
 set output 'test.svg'
 set terminal postscript eps \
 size    5.00000000      \
 cm,    5.00000000     cm \
 enhanced color \
 font 'Helvetica \
 ,           12  \
 ' \
 linewidth    3.00000000    
 set output 'test.eps'
