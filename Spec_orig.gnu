set terminal postscript eps enhanced color 24
        unset label 
        set style line 1 lt -1 lw 1.5 lc rgbcolor "red" pt 1 ps 1 
        set style line 2 lt -1 lw 1.5 lc rgbcolor "blue" pt 2 ps 1 
        set style line 3 lt -1 lw 1.5 lc rgbcolor "web-green" pt 3 ps 1 
        set style line 4 lt -1 lw 1.5 lc rgbcolor "black" pt 4 ps 1 
        set style line 5 lt -1 lw 1.5 lc rgbcolor "magenta" pt 5 ps 1 
        set style line 6 lt -1 lw 1.5 lc rgbcolor "cyan" pt 6 ps 1 
        set style line 7 lt -1 lw 1.5 lc rgbcolor "purple" pt 7 ps 1 
        set style line 8 lt -1 lw 1.5 lc rgbcolor "orange" pt 8 ps 1 
        set style line 9 lt -1 lw 1.5 lc rgbcolor "light-blue" pt 9 ps 1 
        set style line 10 lt -1 lw 1.5 lc rgbcolor "dark-khaki" pt 10 ps 1 
        set style line 11 lt -1 lw 1.5 lc rgbcolor "grey" pt 11 ps 1 
        set style line 12 lt -1 lw 1.5 lc rgbcolor "gold" pt 12 ps 1 
        set style line 13 lt -1 lw 1.5 lc rgbcolor "coral" pt 13 ps 1
        #set label "(b)" right at graph 0.1, graph 0.95
        set xtics 
        set ytics
        #set yrange [*:-8]
        set log y
        set xlabel "s" 
        set ylabel "Energy" 
        set output "test_example.eps"
        set key right bottom
        #plot for [i=2:5] "test_spec.txt" u 1:i w l ls i-1 lw 2 notitle, "test10.txt" u 1:3 w p ls 5 pt 7 title "T_A=10","test100.txt" u 1:3 every 20 w p ls 6 pt 7 title "T_A=100","test1000.txt" u 1:3 every 100 w p ls 7 pt 7 title "T_A=1000",
        plot "test_spec.txt" u 1:($3-$2) w l notitle
