 set terminal pngcairo \
 size          500 ,          500  \
 enhanced font \
 'Helvetica \
 ,           12  \
 '
 set output 'test.png'
 n = 0
 do for [i =            0 :          12 ] {
 n = n + 1
 set output sprintf('tmp/test%d.png',n)
 plot 'test.dat' every ::           0 ::i w l, \
 'test.dat' every ::i::i w p
 }
