# Simulating the two-dimensional time-dependent Schrödinger equation (Project 5)

Our results are summarized in a [project report](report_project1.pdf)

## Desctiption
This folder is contains source code, datafiles and figures for project 5 in FYS4150/3150. This project aims to study the behaviour of a quantum particle in a box by using the normalizes Schrödinger equation. We conduct several experiments using this simulation, for example the single-, double- and triple-slit experiment. Instructions on compiling and running the function, as well as plotting the figures can be found at the bottom of this README file.

![Problem_7 20_outputCube_slits_2 dat_animation](https://user-images.githubusercontent.com/31341364/146273401-9a893a1d-75e6-4bcf-ad0c-5c643b64598e.gif)
![Problem_7 10_outputCube_slits_0 dat_animation](https://user-images.githubusercontent.com/31341364/146273556-28c7a013-3c67-46e6-b9b9-2070487ca593.gif)

## Files
### `main.cpp`
The main file is used to set up and simulate the the different scenarios, by calling on functions from the `src/Crank.cpp` source file.

### `main.exe`
This is the compiled and runnable program. Note that you should run this using the makefile specified at he bottom of this document.

## Folder structure

### [inputs](inputs/)
This folder contains files with input arguments for `main.exe`, and are used to run simulations with different parameters.

### [src](src/)
This folder contains the source code `Crank.cpp` to our program, which is linked during compilation.

### [plotfiles](plotfiles/)
This folder contains the python scrips that are used to plot the animations and figures we have used in the report. Some of these plotfiles also contain some minor calculations performed on the final result. The python scrips are called on when using the makefile and finds the necessary datasets automatically.

### [figures](figures/)
This folder contains all the plots used in the report, which are made by the python scripts. It also contains animations of the different experiments we performed.

### [include](include/)
This folder contains the headerfile `Crank.hpp` for out source file `Crank.hpp`, which is linked in compilation.

### [datafiles](datafiles/)
This folder is the default save location for all the datafiles created when running the program. Most of the datafiles are binary files which are read by the pythonscripts using 'pyarma'. As the files are very large, they are not pushed to git, and you will have to run the program yourself if you want to see them.

## Compiling, running and plotting
The program had some high level parallelisation, but might still take a while to run. The datafiles produced by the program are also pretty large, and are therefore not pushed to the repo. To be able to plot you must then first run the program. Below are instructions on compiling the program, note that you must be in the `project5/` folder for the makefile to run.

### Instructions:
To compile the program, use the following command:

`$ make`

To run the program, use the following command:

`$ make run `

To generate animations and plots, use the following line:

`$ make plot`

To run extra long simulations as specified in `inputs/input.txt`:

`$ make extra`
