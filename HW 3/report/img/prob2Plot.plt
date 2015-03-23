#Setting up the terminal
set terminal postscript enhanced eps color;
#set terminal postscript eps enhanced color \
#    font 'Helvetica,20' linewidth 2
set output 'Prob2Plot.eps';
set size 1,1

# Setting up the axes
set title "Stress-strain behaviour under equibiaxial planar strain"
set xlabel 'F_{11}'
set ylabel 'Components of first Piola-Kirchoff Stress (N)'
set format y "%2.0t{/Symbol \264}10^{%L}"

#For generating legend
set key

#Now plotting
plot "prob2PlotData.dat" using 1:2 title 'P_{11}' with lines lw 6,\
	"prob2PlotData.dat" using 1:3 title 'P_{22}' with lines lw 6;


