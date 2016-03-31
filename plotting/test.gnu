 set terminal pngcairo \
 size          500 ,          500  \
 enhanced font \
 'Helvetica \
 ,           12  \
 '
 set output 'test.png'
 n = 0
 do for [i =            1 :         100 :           1 ] {
 n = n + 1
 set output sprintf('tmp/solar_system%d.png',n)
  splot 'solar_system.dat' \
 u            1 :            2 :            3  \
 every ::           1 ::i w l ls            1  \
 notitle \
 , 'solar_system.dat' u            1 :            2 :            3  \
 every ::i::i w p ls            1  \
 title '' \
 , 'solar_system.dat' \
 u            4 :            5 :            6  \
 every ::           1 ::i w l ls            2  \
 notitle \
 , 'solar_system.dat' u            4 :            5 :            6  \
 every ::i::i w p ls            2  \
 title '' \
 }
