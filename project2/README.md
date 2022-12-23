# Numerical solution of an eigenvalue problem exemplified by a special case of the one-dimensional buckling beam (Project 2).

Our results are summarized in a [project report](report_project2.pdf)

## Description
The main topics of this project was:
- scaling of equations
- eigenvalue problems
- unit testing

We considered a horizontal beam of length L fastened with pin endpoints. The vertical displacement were denoted as u(x), and we were interested in studying the different static shapes the beam (u(x)) can take due to an applied force at the endpoint directed into the beam.

The situation is described as

<img width="183" alt="Screenshot 2022-01-03 at 19 47 13" src="https://user-images.githubusercontent.com/30042718/147967690-8f7f16c5-f877-46e3-a3ef-b893c569469a.png">

which were scaled with x^ = x/L  and &lambda; = (FL<sup>2</sup>)/&gamma;

<img width="153" alt="Screenshot 2022-01-03 at 20 00 14" src="https://user-images.githubusercontent.com/30042718/147968953-ccd80562-8f5a-4eff-ae18-042d814751ea.png">

Declaring v<sub>i</sub> as the discretized approximation to u(x<sub>i</sub>), and inserting boundary conditions v<sub>0</sub> = v<sub>n</sub> = 0, we get the linear algebra eigenvalue problem

<img width="77" alt="Screenshot 2022-01-03 at 20 03 49" src="https://user-images.githubusercontent.com/30042718/147969276-d45a1010-978e-4530-836c-f5680f7a7f4c.png">

where A = tridiag(-1/h<sup>2</sup>, 2/h<sup>2</sup>, -1/h<sup>2</sup>).

The eigenpairs <img width="58" alt="Screenshot 2022-01-03 at 20 07 50" src="https://user-images.githubusercontent.com/30042718/147969692-3e071996-94af-4a9f-94e9-6c4817bd5dae.png"> that solves this equation are then the discretized approximations to the true(analytical) eigenfunctions u<sup>(i)</sup>(x) given by

<img width="474" alt="Screenshot 2022-01-03 at 20 08 48" src="https://user-images.githubusercontent.com/30042718/147969780-7d86bd1e-a939-489b-9efc-237f242630c3.png">

## Compiling, running and plotting

To compile and link main.cpp:

`$ make`

Run the code, main.exe:

`$ make run`

When main.exe is done running you can make plots using the command:

`$ make plot`

Note that the file may take some time to run, as it solves all subtasks as well.
Running main.exe will write results from task 3, 4 and 5 to the terminal. It will also produce some
.txt files with the results from task 6 and 7 stored in the folder *datafiles*

## Folder Structure
Below, you will find a description of each folder. At the bottom you will find instructions on how to compile the program, run it and plot the datafiles that are produced.

### [datafiles](datafiles/)
  This folder contains all the different datafiles(.txt files) that are produced when running the simulation. If the folder is empty, make sure you have compiled and run the program.

### [figures](figures/)
  This folder will contain all the plots that are present in the report. If the folder is empty, make sure you ran the plot command in the terminal.

### [plotfiles](plotfiles/)
  This folder contains the python scrips that were used to plot the results. All the python files saves the figures in the *figure/* folder. Note that if you try to run the program manually, it won't find the correct path to the datafiles, so make sure you run them using the `$make plot` command described above.
