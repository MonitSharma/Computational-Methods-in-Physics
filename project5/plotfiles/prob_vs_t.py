# Problem 7
import numpy as np
import matplotlib.pyplot as plt
import pyarma as pa
import sys

#Adjusting text size
SMALL_SIZE = 13
MEDIUM_SIZE = 14
plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

filename = sys.argv[1]
slits = sys.argv[2]

U_cube = pa.cx_cube() # Creates the complex cube object for storage
U_cube.load("./datafiles/"+str(filename)) #Load the content of the cube you saved into your Python program.

timesteps = U_cube.n_slices

prob_sum = np.array([])
for i in range(0, timesteps): # sum probability for all timesteps
    u = U_cube[pa.single_slice, i] # obtain slice for timestep i
    p = pa.conj(u)@u  # elementwise multiplication to get probability
    psum = np.real(pa.accu(p)) # sum up all probabilites in matrix slice
    prob_sum = np.append(prob_sum, abs(psum-1.0))


t = np.linspace(0,0.008, timesteps)

plt.scatter(t, prob_sum, marker='.')
plt.xlabel("time")
plt.ylabel("|p(x,y) - 1.0|")
plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
plt.yticks([0.0, 0.5*1e-14, 1e-14, 1.5*1e-14, 2.0*1e-14])
plt.ylim([-0.07e-14, 2.3e-14])
plt.grid()
plt.savefig("./figures/prob_vs_t_slits_%s.pdf" %slits)
