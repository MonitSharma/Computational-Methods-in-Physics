# Simulation of the two-dimensional Ising Model (Project 4)

Our results are summarized in a [project report](report_project4.pdf)

## Description
The model will be used to explore temperature-dependent behaviour in ferromagnets and numerically estimate the critical temperature at which our 2D system undergoes a phase transition from a magnetized phase to a phase with no net magnetization.

We estimated the expectancy value of the energy per spin, **< epsilon >**, and of the absolute value of the magnetization per spin, **<|m|>**. This were done by sampling **s** using the Monte Carlo Markov Chain (MCMC) method and taking the mean of **epsilon** and **|m|** over MCMC cycles.
These were the used to find estimations for the specific heat capacity, **C__v**, and the susceptibility, **chi**. These were both normalized to per spin.

We parallelized our code using OpenMP to decrease runtime for large state configuration matrices.  


For full Project description, visit https://anderkve.github.io/FYS3150/book/projects/project4.html

## Folder Structure
Below, you will find a description of each folder. At the bottom you will find instructions on how to compile the program, run it and plot the datafiles that are produced.

### [datafiles](datafiles/)
  This folder contains all the different datafiles(.txt files) that are produced when running the simulation. If the folder is empty, make sure you have compiled and run the program.

### [figures](figures/)
  This folder will contain all the plots that are present in the report. If the folder is empty, make sure you ran the plot command in the terminal.

### [include](include/)
  This folder contains the header files for our program

### [plot_files](plot_files/)
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the *figure/* folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the `$make plot` command described below.

### [src](src/)
  This folder contains the `Ising.cpp` source file. Note that the main.cpp source file will be linked with the `Ising.cpp` source file.

## Code structure
  Our code is mainly split into two parts: `main.cpp` and `Ising.cpp`. We mainly use `main.cpp` to run a simulation for each task using the 'simulate()' function in `main.cpp`. 'simulate()' will then use the functions located in `Ising.cpp`. If you wish to run simulation with some specific parameters, you can do so using the last command given in the instructions below.

## Compiling, running and plotting

  Instructions on compiling running and plotting can be found below. If you are on a different operating system than macOS, you might have to edit line 12 in the makefile to say '-fopenmp' instead of 'lomp'. The program takes a while to run, we recommend you grab a coffee while waiting. The data for problem 6 is not in the repo, if you want to create it you have to run the program. The plotting for task6 is also commented out in the `makefile` and you need to remove the comment if you want to plot it. The plot for task6 is however present in the *figure/* folder.

To build the code:  
`$ make`

To run the code:  
`$ make run`

To make plots:  
`$ make plot`

To test the code using chosen parameters:  
`$ ./main.exe <temperature (float)> <lattice side size (integer)> <MCMC cycles (integer)> <unordered lattice: use 0, ordered lattice: use -1 or 1> <output_file_name>`
