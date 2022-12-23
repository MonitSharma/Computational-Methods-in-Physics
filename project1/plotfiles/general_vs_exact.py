import matplotlib.pyplot as plt
import numpy as np
import math

# Plotting for multiple n values
k = 3 # max exponent of 10 number of points
colors = ['r', 'g', 'y', 'b', 'c', 'm', 'r']

for i in range(1, k+1):
    n = 10**i
    exact_data = np.loadtxt('./datafiles/exact_data%d.txt' %n)
    x = np.array(exact_data[:,0])
    u = np.array(exact_data[:,1])
    genapprox_data = np.loadtxt('./datafiles/approx_general%d.txt' %n)
    v = np.array(genapprox_data[:,1])

    vstar = np.array([0]) # boundary point u_0 = 0
    vstar = np.append(vstar, v)
    vstar = np.append(vstar, 0) # appending boundary point u_1 = 0
    plt.plot(x, vstar, ':', color=colors[i-1], label='v^*(x), N=%d' %n)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.xlabel('x')
plt.title('General algorithm approximation vs exact of Poisson eq.')
plt.plot(x, u, color='k', label='Exact u(x)') # Plotting u(x) with highest n
plt.legend()
plt.savefig('./figures/general_vs_exact.pdf')
