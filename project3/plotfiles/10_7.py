import numpy as np
import matplotlib.pyplot as plt
import sys


SMALL_SIZE = 13
MEDIUM_SIZE = 13
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

runs = ""

if len(sys.argv) > 1 and sys.argv[1] == "average":
	runs = "100avg"


data0 = np.loadtxt('./datafiles/problem10_%sf0.1narrow without interactions.txt'%runs, skiprows=1)
data1 = np.loadtxt('./datafiles/problem10_%sf0.1narrow with interactions.txt'%runs, skiprows=1) #100runs


w0 = np.array(data0[:,0])
f0 = np.array(data0[:,1])
plt.plot(w0, f0, label ="w/o")


w1 = np.array(data1[:,0])
f1 = np.array(data1[:,1])
plt.plot(w1, f1, label ="w/")


plt.title("Amplitude f = 0.1")
plt.xlabel("$\omega_V$, [$\omega_V$] = MHz")
plt.ylabel("Fraction of remaining particles")

#plt.legend()
plt.grid()
plt.subplots_adjust(
	top=0.91,
	bottom=0.14,
	left=0.14,
	right=0.955,
	hspace=0.2,
	wspace=0.2
)

plt.savefig('./figures/fract_vs_angfreq_%s_zoom.pdf'%runs)
#plt.show()
