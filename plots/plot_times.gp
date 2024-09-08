# plot_times.gp

# Input file path
datafile = "times_rhs.txt"

# Set output format and terminal type
set terminal pngcairo size 800,600 enhanced font 'Arial, 12'

# Plot titles, axis labels, and legends
set output 'times_rhs_heat.png'
set title "Time for RHS computation by refinement Level"
set xlabel "Refinement"
set ylabel "Time [s]"
set grid
# set logscale y

# Plot 1_core, 2_cores, and 4_cores with respect to refinement
plot datafile using 1:2 with linespoints lt 1 pt 7 lc rgb "red" title "1 Core", \
     datafile using 1:3 with linespoints lt 2 pt 7 lc rgb "blue" title "2 Cores", \
     datafile using 1:4 with linespoints lt 3 pt 7 lc rgb "green" title "4 Cores"
