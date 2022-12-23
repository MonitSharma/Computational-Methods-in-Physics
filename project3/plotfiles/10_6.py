import numpy as np
import matplotlib.pyplot as plt

#read out w_v and N for each f as w_1,w_2,w_3,N_1,N_2,N_3


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


# f = 0.1
data1 = np.loadtxt('./datafiles/problem10_f0.1broad.txt', skiprows=1)
w1 = np.array(data1[:,0])
f1 = np.array(data1[:,1])
plt.plot(w1, f1, label ="f = 0.1")

# f = 0.4
data2 = np.loadtxt('./datafiles/problem10_f0.4broad.txt', skiprows=1)
w2= np.array(data2[:,0])
f2 = np.array(data2[:,1])
plt.plot(w2,f2, label ="f = 0.4")

# f= 0.7
data3 = np.loadtxt('./datafiles/problem10_f0.7broad.txt', skiprows=1)
w3 = np.array(data3[:,0])
f3 = np.array(data3[:,1])
plt.plot(w3,f3, label ="f = 0.7")

plt.xlabel("$\omega_V$, [$\omega_V$] = MHz")
plt.ylabel("Fraction of remaining particles")
plt.legend()
plt.grid()
plt.subplots_adjust(
	top=0.91,
	bottom=0.14,
	left=0.14,
	right=0.955,
	hspace=0.2,
	wspace=0.2
)

plt.savefig('./figures/fraction_vs_angfreq.pdf')
#plt.show()
