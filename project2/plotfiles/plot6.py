import matplotlib.pyplot as plt
from math import *
import numpy as np
import sys


N =[]
itterations =[]

#reads from datasets and inserts into lists
with open("./datafiles/task6_dataset.txt") as myfile: #opens file
    myfile.readline()   #skips first line
    for line in myfile:
        mydata = line.split(",")
        N.append(float(mydata[0]))
        itterations.append(float(mydata[1]))

myfile.close()

#Converts to array
N = np.array(N)
itterations = np.array(itterations)


#plots function
plt.plot(N,itterations, label ="Number of rotations")

#Sets up, saves and shows figure
plt.xlabel("N values")
plt.ylabel("Number of rotations")
plt.title("Number of rotations when A is a symmetrical tridiagonal NxN matrix")
plt.legend()
plt.grid()
plt.savefig('./figures/rotations_vs_N.pdf')
plt.show()
