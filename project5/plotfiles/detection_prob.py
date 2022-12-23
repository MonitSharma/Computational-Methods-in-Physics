# Problem 9
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

#Input arguments that determines what is being plotted.
filename = sys.argv[1]
slits = sys.argv[2]

U = pa.cx_mat() # Create pa.mat object (just as arma::mat in C++)
U.load("./datafiles/"+str(filename)) # Load the content of the matrix you saved into your Python program.
U = np.array(U) # U-matrix at time t=0.002
prob_matrix = np.conj(U)*U # probability matrix


h = 0.005 # stepsize for position
index = int(0.8/h) # index at which x = 0.8

p = np.real(prob_matrix[index]) # get array of probabilities at x = 0.8
p_norm = p/np.sqrt(np.sum(p**2)) # normalise

points = int((1/h) +1) # M = 201
y = np.linspace(0,1,points)

plt.plot(y, p_norm)
plt.xlabel("y [Units of distance/1]")
plt.ylabel("p(y|x=0.8; t=0.002)")
plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
plt.yticks([0.0, 0.1, 0.2, 0.3])
plt.ylim([-0.01, 0.3])
plt.grid()
plt.savefig("./figures/detection_prob_slits_%s.pdf" %slits)
