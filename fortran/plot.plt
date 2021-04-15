
set terminal png truecolor size 1080,1080    
set output "spin.png"    
set autoscale   
set linetype 1 linecolor "blue"
set xlabel "time(ps)"   
set ylabel "Magnetization"    
set title "Spin dynamics"    
set grid    
set key top left box
set key samplen 2
set key width 3
set key font ',13'
set border 3
set border lw 5
set tics nomirror
plot "data1.dat" using 1:2 title "Mx" with linespoints,'' using 1:3 title "My" with linespoints,'' using 1:4 title "Mz" with linespoints
