 set terminal pngcairo \
 size          800 ,          800  \
 enhanced font \
 'Helvetica \
 ,            1  \
 '
 n = 0
 do for [i =            0 :         345 :           1 ] {
 n = n + 1
 set output sprintf('tmp/solar_system%d.png',n)
  splot 'solar_system.dat' \
 u            1 :            2 :            3  \
 every ::           0 ::i w l ls            1  \
 , 'solar_system.dat' u            1 :            2 :            3  \
 every ::i::i w p ls            1  \
 , 'solar_system.dat' \
 u            4 :            5 :            6  \
 every ::           0 ::i w l ls            2  \
 , 'solar_system.dat' u            4 :            5 :            6  \
 every ::i::i w p ls            2  \
 , 'solar_system.dat' \
 u            7 :            8 :            9  \
 every ::           0 ::i w l ls            3  \
 , 'solar_system.dat' u            7 :            8 :            9  \
 every ::i::i w p ls            3  \
 , 'solar_system.dat' \
 u           10 :           11 :           12  \
 every ::           0 ::i w l ls            4  \
 , 'solar_system.dat' u           10 :           11 :           12  \
 every ::i::i w p ls            4  \
 , 'solar_system.dat' \
 u           13 :           14 :           15  \
 every ::           0 ::i w l ls            5  \
 , 'solar_system.dat' u           13 :           14 :           15  \
 every ::i::i w p ls            5  \
 }
