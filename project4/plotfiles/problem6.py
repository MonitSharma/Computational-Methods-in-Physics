# Estimating the pdf of epsilon
import numpy as np
#from plotlib import makeplots
import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

#Adjusting text size
SMALL_SIZE = 13
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

fig, ax = plt.subplots(figsize = (7, 6))
plt.subplots_adjust(
top=0.95,
bottom=0.15,
left=0.15,
right=0.95,
hspace=0.2,
wspace=0.2
)

for i in range(0,2):
    if i == 0:
        T = 1.0
        alpha = 1 # transparency factor
        bins = 10
    else:
        T = 2.4
        alpha = 0.5
        bins = 50

    # comment out to experiment with hist look
    data = np.loadtxt('./datafiles/histogram_T_%.1f_unordered.txt' %T, skiprows=1)
    eps = np.array(data[:])

    print("\n---------- T = %.1f ----------" %T)
    print("n_samples_total = ", len(eps))

    eps_mean = np.mean(eps)
    print("eps_mean = ", eps_mean)
    eps_variance = np.sum((eps - eps_mean)**2) / (len(eps)-1)
    print("Var(eps) = ", eps_variance)
    print("------------------------------\n")

    ax.hist(eps, label="T=%.1f"%T , density=True, log=True, alpha=alpha, stacked=True, bins=bins)
    ax.plot([eps_mean, eps_mean], [0.0, ax.get_ylim()[1]], linestyle="dashed", linewidth=1.5, label="sample mean")

# Setting ticks on x-axis
ax.xaxis.set_major_locator(MultipleLocator(0.20))
ax.xaxis.set_major_formatter('{x:.2f}')
ax.xaxis.set_minor_locator(MultipleLocator(0.01))

plt.ylabel('$p_{\epsilon}(\epsilon)_{est}$')
plt.xlabel('$\epsilon$')
plt.grid()
plt.legend(loc="upper right")
plt.savefig("./figures/histogram.pdf")
plt.show()
