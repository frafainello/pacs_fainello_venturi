# plots.gp

# Read the input filename from the command-line argument
filename = ARG1

# Check if the filename was provided
if (strlen(filename) == 0) {
    print "Usage: gnuplot -c plots.gp 'your_data_file.txt'"
    exit
}

# Set output file type and font
set terminal pngcairo size 800,600 enhanced font 'Arial, 12'

# Plot 1: Iterations vs residual
set output 'plots/plot_max_min_z.png'
set title "Iterations vs residual"
set xlabel "Iterations"
set ylabel "Residual"
set grid
plot filename using 1:2 with linespoints title 'max', \
     filename using 1:3 with linespoints title 'min'

# Plot 2: Iterations vs gradient norms
set output 'plots/plot_max_min_grad.png'
set title "Iterations vs gradient norms"
set xlabel "Iterations"
set ylabel "Gradient Norms"
set grid
plot filename using 1:4 with linespoints title 'max', \
     filename using 1:5 with linespoints title 'min'
