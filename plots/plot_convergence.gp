# plot_convergence.gp

# Input file paths
file1 = "vtk_files/mesh_10/metrics/np4_face_heat_standard.txt"
file2 = "vtk_files/mesh_10/metrics/np4_face_heat_penalty.txt"
file3 = "vtk_files/mesh_10/metrics/np4_face_heat_lagrangian.txt"

title1 = 'Standard'
title2 = 'Penalty'
title3 = 'ADMM'

# Set output format and terminal type
set terminal pngcairo size 800,600 enhanced font 'Arial, 12'


# Plot titles, axis labels, and legends
set output 'plots/convergence_z.png'
set title "Convergence plot for increment"
set xlabel "Iterations"
set ylabel "Increment"
set yrange [-0.003:0.003]
set grid

# First plot: iterations vs max_z and min_z
plot file1 using 1:2 with linespoints lt 1 pt 7 lc rgb "red" title sprintf("%s - max", title1), \
     file1 using 1:3 with linespoints lt 1 pt 5 lc rgb "red" title sprintf("%s - min", title1), \
     file2 using 1:2 with linespoints lt 2 pt 7 lc rgb "blue" title sprintf("%s - max", title2), \
     file2 using 1:3 with linespoints lt 2 pt 5 lc rgb "blue" title sprintf("%s - min", title2), \
     file3 using 1:2 with linespoints lt 3 pt 7 lc rgb "green" title sprintf("%s - max", title3), \
     file3 using 1:3 with linespoints lt 3 pt 5 lc rgb "green" title sprintf("%s - min", title3)


# Second plot: iterations vs max_grad and min_grad
set output 'plots/convergence_grad.png'
set title "Convergence plot for gradient norms"
set xlabel "Iterations"
set ylabel "Gradient norms"
set yrange [0.8:1.2]
set grid

plot file1 using 1:4 with linespoints lt 1 pt 7 lc rgb "red" title sprintf("%s - max", title1), \
     file1 using 1:5 with linespoints lt 1 pt 5 lc rgb "red" title sprintf("%s - min", title1), \
     file2 using 1:4 with linespoints lt 2 pt 7 lc rgb "blue" title sprintf("%s - max", title2), \
     file2 using 1:5 with linespoints lt 2 pt 5 lc rgb "blue" title sprintf("%s - min", title2), \
     file3 using 1:4 with linespoints lt 3 pt 7 lc rgb "green" title sprintf("%s - max", title3), \
     file3 using 1:5 with linespoints lt 3 pt 5 lc rgb "green" title sprintf("%s - min", title3)
     
set output
