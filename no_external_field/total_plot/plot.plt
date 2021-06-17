
set terminal png truecolor size 1080,1080    
set output "spin.png"    
set autoscale   
set linetype 1 linecolor "blue" lw 0.01
set xlabel "time(ps)"   
set ylabel "Magnetization"    
set title "Spin_dynamics"    
set grid    
set key top right
set key samplen 2
set key width 3
#set key font ',20'
set border 
set border lw 5
set tics nomirror
set yrange [0.4:1.2]
plot "moment_tot.dat" using 1:4 title "Mx" with linespoints,'' using 1:5 title "My" with linespoints,'' using 1:6 title "Mz" with lines lw 3
