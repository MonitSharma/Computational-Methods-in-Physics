# Simulating the Behaviour of Charged Particles in a Penning trap (Project 3)
Our results are summarized in the [project report](report_project3.pdf).

![Penning trap simulation in spacial plane](https://github.com/cecilmkl/comp_phys/blob/main/project3/figures/xy_1000.png)

## Description
A Penning trap is a device that utilizes a static configuration of electic and magnetic fields to confine charged particles.

## Build, run and plot
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
  This folder will contain all the plots that are present in the report created from the Python files in the *plot_files/* folder. If the folder is empty, make sure to run the plot command in the terminal.

### [plot_files](plot_files/)
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the *figure/* folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the `$ make plot` command described above.

### [include](include/)
  This folder contains the header files for our program; `Particle.hpp` and `PenningTrap.hpp`

### [src](src/)
  This folder contains the `Particle.cpp` and `PenningTrap.cpp` source file. Note that the `main.cpp` source file will be linked with these files.
