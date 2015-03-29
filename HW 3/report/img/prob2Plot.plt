#Setting up the terminal
set terminal postscript enhanced eps color;
#set terminal postscript eps enhanced color \
#    font 'Helvetica,20' linewidth 2
set output 'EquiBiaxial.eps';
set size 1,1

# Setting up the axes
set title "Stress-strain Behaviour under Equibiaxial Planar Strain"
set xlabel 'F_{11}'
set ylabel 'Components of first Piola-Kirchoff Stress (N)'
set format y "%2.0t{/Symbol \264}10^{%L}"

#For generating legend
set key

#Now plotting
plot "EquibiaxialStrainData.dat" using 1:2 title 'P_{11}' with lines lw 6,\
	"EquibiaxialStrainData.dat" using 1:3 title 'P_{22}' with lines lw 6;


