 set terminal wxt size 10,2 enhanced font 'Helvetica,14' persist
 set terminal pngcairo size 10,2 enhanced font 'Helvetica,14'
 set output 'test.png'
 set terminal svg size 10,2 fname 'Helvetica fsize 14
 set output 'test.svg'
 set terminal postscript eps size 10,2 enhanced color font 'Helvetica,14' linewidth 2
 set output 'test.eps'
 set terminal epslatex size 10,2 color colortext 14
 set output 'test.tex'
 plot 'test.dat' u   1:   2, 'test.dat' u   3:   4, 'test.dat' u   5:   6
