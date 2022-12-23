## Folder content
This folder contains the python scipts used for plotting. Alle the plotfiles save their figure to the `project5/figures/` folder.

## Scripts
### `Animation.py`
This scrips creates the animation of the simulation of both the probability function, and the real and imaginary part of u. It also creates snapshots of the animation if requested in the input arguments. These are specified in the `makefile` in the `project5/` folder. Note what this simulation aslo calculates the complex conjugate of u, as we had some trouble with getting compex conjugate in cpp to work.  

### `detection_prob.py`
This script creates the one dimensional plot of the probability of detecting the particle along a virtual detection screen at x=0.8. It takes the filename of the datafile as well as the number of slits the data was produced for as input arguments. The output is a pdf of the plot saved in the `project5/figures/` folder.

### `prob_vs_t.py`
This script creates the one dimensional plot of the deviation from 1.0 for the total probability as a function of time. It takes the filename of the datafile as well as the number of slits the data was produced for as input arguments. The output is a pdf of the plot saved in the `project5/figures/` folder.

