 set terminal pngcairo \
 size          800 ,          800  \
 enhanced font \
 'Helvetica \
 ,           12  \
 '
 set xrange [  -1.65026855     :    1.39538169     ]
 set yrange [  -1.45760167     :    1.58021462     ]
 set zrange [  -5.12485579E-02 :    5.33940010E-02 ]
 set xlabel 'X, A.U.'
 set ylabel 'Y, A.U.'
 set zlabel 'Z, A.U.'
 set title 'Inner Solar System'
 set xtics   0.500000000    
 set mztics            2
 n = 0
 do for [i =            0 :         346 :           1 ] {
 n = n + 1
 set output sprintf('tmp/solar_system%d.png',n)
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
