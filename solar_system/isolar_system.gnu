 set terminal pngcairo \
 size         1000 ,         1000  \
 enhanced font \
 'TeX Gyre Pagella \
 ,           13  \
 '
 set xrange [  -1.79999995     :    1.79999995     ]
 set yrange [  -1.79999995     :    1.79999995     ]
 set zrange [  -5.99999987E-02 :    5.99999987E-02 ]
 set xlabel 'x, A.U.'
 set ylabel 'y, A.U.'
 set zlabel 'z, A.U.'
 set title 'Inner Solar System'
 set xtics   -1.79999995     ,   0.400000006     ,    1.79999995    
 set ytics   -1.79999995     ,   0.400000006     ,    1.79999995    
 set ztics   -5.99999987E-02 ,    1.99999996E-02 ,    5.99999987E-02
 set mxtics            2
 set mytics            2
 set mztics            2
 set grid xtics
 set grid ytics
 set grid ztics
 set xyplane at  -5.99999987E-02
 n = 0
 do for [i =            0 :         346 :           1 ] {
 n = n + 1
 set output sprintf('tmp/isolar_system%d.png',n)
 splot 'solar_system.dat' \
 u            1 :            2 :            3  \
 every ::           0 ::i w l ls            1  \
 notitle \
 , 'solar_system.dat' u            1 :            2 :            3  \
 every ::i::i w p ls            1  \
 title 'Sun' \
 , 'solar_system.dat' \
 u            4 :            5 :            6  \
 every ::           0 ::i w l ls            2  \
 notitle \
 , 'solar_system.dat' u            4 :            5 :            6  \
 every ::i::i w p ls            2  \
 title 'Mercury' \
 , 'solar_system.dat' \
 u            7 :            8 :            9  \
 every ::           0 ::i w l ls            3  \
 notitle \
 , 'solar_system.dat' u            7 :            8 :            9  \
 every ::i::i w p ls            3  \
 title 'Venus' \
 , 'solar_system.dat' \
 u           10 :           11 :           12  \
 every ::           0 ::i w l ls            4  \
 notitle \
 , 'solar_system.dat' u           10 :           11 :           12  \
 every ::i::i w p ls            4  \
 title 'Earth' \
 , 'solar_system.dat' \
 u           13 :           14 :           15  \
 every ::           0 ::i w l ls            5  \
 notitle \
 , 'solar_system.dat' u           13 :           14 :           15  \
 every ::i::i w p ls            5  \
 title 'Mars' \
 }
