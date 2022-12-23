
import numpy as np
#from plotlib import makeplots
import matplotlib.pyplot as plt


#Adjusting text size
SMALL_SIZE = 14
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


data_units =  ['$N_{cycles}$','E','M', r'$[<\varepsilon>]$=J', '[<|m|>] = 1', r'$[C_v] = k_B^2$', r"[$\chi$] = $\frac{1}{k_B}$"]
data_labels = ['MC cycles','E','M', 'eps', 'm', 'C_v', "susceptibility"]
data = np.loadtxt('./datafiles/task4.txt',skiprows=1)
analytical_values = np.loadtxt('./datafiles/analytical_2x2_T=1.0.txt',skiprows=1)


for i in range(3,7):
	fig, ax = plt.subplots(figsize = (6, 5))
	plt.subplots_adjust(
	top=0.95,
	bottom=0.15,
	left=0.2,
	right=0.95,
	hspace=0.2,
	wspace=0.2
	)
	plt.plot(np.array(data[:,0]),np.array(data[:,i]), label="MCMC")
	plt.plot(np.array(data[:,0]),np.ones(len(data[:,0]))*analytical_values[i-2], label='Analytical')
	#ax.set_title("Training mse")
	ax.set_ylabel(data_units[i])
	ax.set_xlabel(data_units[0])
	plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
	plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
	plt.grid()
	plot_name = "Task4_"+data_labels[i]+"_figure.pdf"
	plt.legend()
	plt.savefig('./figures/'+plot_name)
	plt.clf()
	plt.cla()
