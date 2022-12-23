import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Plot of z-trajectory of 1 particle

data = np.loadtxt('./datafiles/RK4_i_10000_d_100_p_1_pi_1_outputs_tzv_pert_0_rs_0_f_0.0_w_v_0.0.txt', skiprows=1)
t = np.array(data[:,0])
z = np.array(data[:,1])



SMALL_SIZE = 17
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# Tick frequency on axes
ax = plt.axes()
ax.xaxis.set_major_locator(ticker.MultipleLocator(10)) # 5 = tick increment


plt.plot(t, z)
plt.xlabel('time, [t] = $\mu s$')
plt.ylabel('Position z(t), [z] = $\mu m$')
plt.grid()

plt.subplots_adjust(
	top=0.925,
	bottom=0.135,
	left=0.195,
	right=0.98,
	hspace=0.2,
	wspace=0.2
)


plt.savefig('./figures/zt.pdf')
#plt.show()
