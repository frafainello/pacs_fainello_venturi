# plot_iterations.gp

# Input file path
file = "vtk_files/iterations.txt"

title1 = 'Standard'
title2 = 'Penalty'
title3 = 'Lagrangian'

# Set output format and terminal type
set terminal pngcairo size 800,600 enhanced font 'Arial, 12'

# Plot titles, axis labels, and legends
set output 'plots/iteration.png'
set title "Iterations for convergence"
set xlabel "Refinement"
set ylabel "Iterations"
set grid

# Plot data from the file
plot file using 1:2 with linespoints lt 1 pt 7 lc rgb "red" title title1, \
     file using 1:3 with linespoints lt 2 pt 7 lc rgb "blue" title title2, \
     file using 1:4 with linespoints lt 3 pt 7 lc rgb "green" title title3

set output
