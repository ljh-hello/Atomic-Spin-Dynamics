#!/bin/bash
# for jj in {0.01,0.03,0.05,0.07,0.1,0.2,0.3,0.4,0.5}
for jj in {0.01,0.5}
do
mkdir damping_size_$jj
cp damping_test.f90 plot.plt damping_size_$jj
cd damping_size_$jj
sed -i '' "s/damping_onsite=0.1/damping_onsite=$jj/g" damping_test.f90
mkdir total_plot
cp plot.plt total_plot
cd total_plot
sed -i '' '/moment_tot.dat/d' ./plot.plt
sed -i  '' "s/Spin_dynamics/Spin dynamics/g "  plot.plt
cd ..
for ii in {1..5}
do
mkdir damp_type_$ii
cp damping_test.f90 plot.plt damp_type_$ii
cd damp_type_$ii

sed -i '' "s/d_type=1/d_type=$ii/g" damping_test.f90 
gfortran damping_test.f90 
./a.out
if [ $ii == 1 ]
then 
sed -i  '' "s/Spin_dynamics/Diagonal onsite damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Diagonal-onsite-damping.dat
# echo "plot \"../damp_type_$ii/moment_tot.dat\" using 1:6 title \"Diagonal onsite damping\" with linespoints" >> ../total_plot/plot.plt
elif [ $ii == 2 ]
then 
sed -i  '' "s/Spin_dynamics/Different onsite damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Different-onsite-damping.dat
elif [ $ii == 3 ]
then 
sed -i  '' "s/Spin_dynamics/Full onsite damping/g  "  plot.plt
cp moment_tot.dat ../total_plot/Full-onsite-damping.dat

elif [ $ii == 4 ]
then 
sed -i  '' "s/Spin_dynamics/Diagonal nonlocal damping/g  "  plot.plt
cp moment_tot.dat ../total_plot/Diagonal-nonlocal-damping.dat

elif [ $ii == 5 ]
then 
sed -i '' "s/Spin_dynamics/Full nonlocal damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Full-nonlocal-damping.dat

fi

# gnuplot plot.plt
printf $jj,$ii \n
cd ..
done
cd total_plot
echo "filenames = \"Diagonal-onsite-damping Different-onsite-damping Full-onsite-damping Diagonal-nonlocal-damping Full-nonlocal-damping\""  >> plot.plt
echo "plot for [file in filenames] file.\".dat\" using 1:6 with lines lw 5" title "file"  >> plot.plt
gnuplot plot.plt
cd ..
cd ..
done
