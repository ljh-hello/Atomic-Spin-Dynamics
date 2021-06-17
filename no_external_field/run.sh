#!/bin/bash
mkdir total_plot
cp plot.plt total_plot
cd total_plot
#sed -i '' '/moment_tot.dat/d' ./plot.plt
#sed -i  '' "s/Spin_dynamics/Spin dynamics with different non-local damping(B_{ext}=0)/g "  plot.plt
cd ..
for ii in 4
do
mkdir damp_type_$ii
cp damping_test.f90 POSCAR2 plot.plt damp_type_$ii
cd damp_type_$ii
mkdir total_plot
sed -i '' "s/d_type=1/d_type=$ii/g" damping_test.f90 
#gfortran damping_test.f90 
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

elif [ $ii == 5 ]
then 
sed -i '' "s/Spin_dynamics/Full nonlocal damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Full-nonlocal-damping.dat

elif [ $ii == 6 ]
then 
sed -i '' "s/Spin_dynamics/Negative full nonlocal damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Negative-full-nonlocal-damping.dat

elif [ $ii == 7 ]
then 
sed -i '' "s/Spin_dynamics/Negative diagonal nonlocal damping/g "  plot.plt
cp moment_tot.dat ../total_plot/Negative-diagonal-nonlocal-damping.dat

fi

for jj in {0.1,0.3,0.5,0.7,1}
do
mkdir nonlocal_scale_$jj
cp damping_test.f90 POSCAR2  nonlocal_scale_$jj
cd nonlocal_scale_$jj
sed -i '' "s/non_scale=0.3/non_scale=$jj/g" damping_test.f90
gfortran damping_test.f90 
./a.out
cp moment_tot.dat ../total_plot/$jj.dat
# gnuplot plot.plt
printf $ii,$jj
cd ..
done

echo "filenames = \"0.1 0.3 0.5 0.7 1 \""  >> total_plot/plot.plt
echo "plot for [file in filenames] file.\".dat\" using 1:6 with lines lw 5" title "file"  >> total_plot/plot.plt
cd total_plot
gnuplot plot.plt
cd ..
done
#cd total_plot
#echo "filenames = \"Diagonal-onsite-damping Different-onsite-damping Full-onsite-damping Diagonal-nonlocal-damping Full-nonlocal-damping Negative-full-nonlocal-damping Negative-diagonal-nonlocal-damping\""  >> plot.plt
#echo "plot for [file in filenames] file.\".dat\" using 1:6 with lines lw 5" title "file"  >> plot.plt
#gnuplot plot.plt
cd ..
