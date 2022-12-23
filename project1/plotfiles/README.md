# Folder for all plotting files

## plot_exact.py

Python script that reads the data in exact_data{N}.txt and plots the exact solution of the Poisson equation, u(x).
Output: 'figures/exact.pdf'

## general_vs_exact.py

Plots the general approximations from approx_general<N>.txt against the exact values exact_data<N>.txt.
Output: 'figures/general_vs_exact.pdf'


## errorplot.py

Plots logarithmic10 absolute and relative errors from file errors<N>.txt against x-values.
Output:'figures/error_plot_log.pdf'
