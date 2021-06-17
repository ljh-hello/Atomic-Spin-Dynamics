
set terminal png truecolor size 1080,1080    
set output "spin.png"    
set autoscale   
set linetype 1 linecolor "blue" lw 0.01
set xlabel "time(ps)"   
set ylabel "Magnetization"    
set title "Diagonal nonlocal damping"    
set grid    
set key top right
set key samplen 2
set key width 3
set key font ',20'
set border 
set border lw 5
set tics nomirror
set yrange [0.5:1.1]
set xrange [0:0.8]
filenames = "0.1 0.3 0.5 0.7 1 "
plot for [file in filenames] file.".dat" using 1:6 with lines lw 5 title file
