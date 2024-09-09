# plot_convergence.gp

# Input file paths
file1 = "../vtk_files/mesh_10/metrics/np4_face_heat_standard.txt"
file2 = "../vtk_files/mesh_10/metrics/np4_face_heat_penalty.txt"
file3 = "../vtk_files/mesh_10/metrics/np4_face_heat_lagrangian.txt"

title1 = 'Standard'
title2 = 'Penalty'
title3 = 'ADMM'

# Set output format and terminal type
set terminal pngcairo size 800,600 enhanced font 'Arial, 12'

# Define line style for dashed black line
set style line 1 lt 2 lc rgb "black" lw 2 dashtype 2 # Dashed line style

# Plot titles, axis labels, and legends
set output 'convergence_z_heat.png'
set title "Convergence plot for increment"
set xlabel "Iterations"
set ylabel "Increment"
set yrange [*:0.1]
set grid

# First plot: iterations vs max_z and min_z
plot file1 using 1:2 with linespoints lt 1 pt 7 lc rgb "red" title title1, \
     file2 using 1:2 with linespoints lt 2 pt 7 lc rgb "blue" title title2, \
     file3 using 1:2 with linespoints lt 3 pt 7 lc rgb "green" title title3, \
     0 with lines ls 1 notitle # Add dashed line at y = 0


# Second plot: iterations vs max_grad and min_grad
set output 'convergence_grad_heat.png'
set title "Convergence plot for gradient norms"
set xlabel "Iterations"
set ylabel "Gradient norms"
set yrange [0.8:1.2]
set grid

plot file1 using 1:3 with linespoints lt 1 pt 7 lc rgb "red" title sprintf("%s - max", title1), \
     file1 using 1:4 with linespoints lt 1 pt 5 lc rgb "red" title sprintf("%s - min", title1), \
     file2 using 1:3 with linespoints lt 2 pt 7 lc rgb "blue" title sprintf("%s - max", title2), \
     file2 using 1:4 with linespoints lt 2 pt 5 lc rgb "blue" title sprintf("%s - min", title2), \
     file3 using 1:3 with linespoints lt 3 pt 7 lc rgb "green" title sprintf("%s - max", title3), \
     file3 using 1:4 with linespoints lt 3 pt 5 lc rgb "green" title sprintf("%s - min", title3), \
     1 with lines ls 1 notitle # Add dashed line at y = 1

set output
