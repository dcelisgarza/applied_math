 set terminal pngcairo \
 size         1000 ,         1000  \
 enhanced font \
 'TeX Gyre Pagella \
 ,           13  \
 '
 set xrange [  -31.0000000     :    44.0000000     ]
 set yrange [  -33.0000000     :    45.0000000     ]
 set zrange [  -15.0000000     :    9.00000000     ]
 set xlabel 'x, A.U.'
 set ylabel 'y, A.U.'
 set zlabel 'z, A.U.'
 set title 'Outer Solar System'
 set xtics   -31.0000000     ,    15.0000000     ,    44.0000000    
 set ytics   -33.0000000     ,    10.0000000     ,    47.0000000    
 set ztics   -15.0000000     ,    4.00000000     ,    9.00000000    
 set mxtics            2
 set mytics            2
 set mztics            2
 set grid xtics
 set grid ytics
 set grid ztics
 set xyplane at  -15.0000000    
 n = 0
 do for [i =            0 :         346 :           5 ] {
 n = n + 1
 set output sprintf('tmp/osolar_system%d.png',n)
 splot 'solar_system.dat' \
 u           16 :           17 :           18  \
 every ::           0 ::i w l ls            1  \
 notitle \
 , 'solar_system.dat' u           16 :           17 :           18  \
 every ::i::i w p ls            1  \
 title 'Jupiter' \
 , 'solar_system.dat' \
 u           19 :           20 :           21  \
 every ::           0 ::i w l ls            2  \
 notitle \
 , 'solar_system.dat' u           19 :           20 :           21  \
 every ::i::i w p ls            2  \
 title 'Saturn' \
 , 'solar_system.dat' \
 u           22 :           23 :           24  \
 every ::           0 ::i w l ls            3  \
 notitle \
 , 'solar_system.dat' u           22 :           23 :           24  \
 every ::i::i w p ls            3  \
 title 'Uranus' \
 , 'solar_system.dat' \
 u           25 :           26 :           27  \
 every ::           0 ::i w l ls            4  \
 notitle \
 , 'solar_system.dat' u           25 :           26 :           27  \
 every ::i::i w p ls            4  \
 title 'Neptune' \
 , 'solar_system.dat' \
 u           28 :           29 :           30  \
 every ::           0 ::i w l ls            5  \
 notitle \
 , 'solar_system.dat' u           28 :           29 :           30  \
 every ::i::i w p ls            5  \
 title 'Pluto' \
 }
