#!/bin/bash

# Define the Gnuplot files to run
gnuplot_files=("plot_convergence.gp" "plot_iterations.gp" "plot_times.gp" "plot_convergence_ic_like_bc.gp")

# Loop through each file and run it with Gnuplot
for file in "${gnuplot_files[@]}"
do
    echo "Running $file..."
    gnuplot "$file"
done

echo "All Gnuplot scripts have been executed."
