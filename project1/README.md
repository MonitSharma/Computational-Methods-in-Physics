# Numerical solution to the one-dimensional Poisson equation by converting the second derivative to a solvable matrix equation Av = g (Project 1)

Our results are summarized in a [project report](report_project1.pdf)

## Description

`main.cpp` computes the exact solution, u(x), and the numerical approximations v(x) and v<sup>*</sup>(x) using both the general algorithm and special algorithm for solving the matrix equation.
This is done with different number of steps, N = 10<sup>i</sup> for i = 1,2,...,6

The code outputs the result to files based on N:

  `exact_data{N}.txt`

  `approx_general{N}.txt`

  `approx_special{N}.txt`

It also stores log10 absolute and relative error, as well as max(rel_error) in files:

  `errors{N}.txt`

The number of steps can be changed manually in the makefile, or by running:

`$ ./main.cpp {N}`

## Compiling, running and plotting

To build the code:  
`$ make`

To run the code:  
`$ make run`

To make plots:  
`$ make plot`


## Folder Structure

Below, you will find a description of each folder. At the bottom you will find instructions on how to compile the program, run it and plot the datafiles that are produced.

### [datafiles](datafiles/)
  This folder contains all the different datafiles(.txt files) that are produced when running the simulation. If the folder is empty, make sure you have compiled and run the program.

### [figures](figures/)
  This folder will contain all the plots that are present in the report. If the folder is empty, make sure you ran the plot command in the terminal.

### [plotfiles](plotfiles/)
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the `figures/`` folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the '$make plot' command described above.
