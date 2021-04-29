#!/bin/bash

for ii in 5
do
mkdir damp_type_$ii
cp damping_test.f90 plot.plt damp_type_$ii
cd damp_type_$ii

sed -i '' "s/d_type=1/d_type=$ii/g" damping_test.f90 
if [ $ii == 1 ]
then 
sed -i  '' "s/Spin_dynamics/Diagonal onsite damping/g "  plot.plt

# echo "plot \"../damp_type_$ii/moment_tot.dat\" using 1:6 title \"Diagonal onsite damping\" with linespoints" >> ../total_plot/plot.plt
elif [ $ii == 2 ]
then 
sed -i  '' "s/Spin_dynamics/Different onsite damping/g "  plot.plt

elif [ $ii == 3 ]
then 
sed -i  '' "s/Spin_dynamics/Full onsite damping/g  "  plot.plt


elif [ $ii == 4 ]
then 
sed -i  '' "s/Spin_dynamics/Diagonal nonlocal damping/g  "  plot.plt


elif [ $ii == 5 ]
then 
sed -i '' "s/Spin_dynamics/Full nonlocal damping/g "  plot.plt

elif [ $ii == 6 ]
then 
sed -i '' "s/Spin_dynamics/Negative full nonlocal damping/g "  plot.plt

sed -i '' "s/step=800000/step=1300000/g "  damping_test.f90 

elif [ $ii == 7 ]
then 
sed -i '' "s/Spin_dynamics/Negative diagonal nonlocal damping/g "  plot.plt
sed -i '' "s/step=800000/step=1300000/g "  damping_test.f90 
fi
mkdir total_plot
cp plot.plt total_plot
cd total_plot
sed -i '' '/moment_tot.dat/d' ./plot.plt
cd ..
for jj in {0.01,0.03,0.05,0.07,0.1,0.3,0.5,0.7,1}
do
mkdir nonlocal_scale_$jj
cp damping_test.f90 nonlocal_scale_$jj
cd nonlocal_scale_$jj
sed -i '' "s/non_scale=0.5/non_scale=$jj/g" damping_test.f90
gfortran damping_test.f90 
./a.out
cp moment_tot.dat ../total_plot/$jj.dat
# gnuplot plot.plt
printf $ii,$jj
cd ..
done
echo "filenames = \"0.01 0.03 0.05 0.07 0.1 0.3 0.5 0.7 1 \""  >> total_plot/plot.plt
echo "plot for [file in filenames] file.\".dat\" using 1:6 with lines lw 5" title "file"  >> total_plot/plot.plt
cd total_plot
gnuplot plot.plt
cd ../..
done