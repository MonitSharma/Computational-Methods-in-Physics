import matplotlib.pyplot as plt
import numpy as np
import math

k = 6 # max exponent of 10 number of points
colors = ['r', 'g', 'y', 'b', 'c', 'm', 'r']

fig, (ax1, ax2)= plt.subplots(2, sharex=True)
fig.suptitle('Errors')

for i in range(1, k+1):
    N = 10**i
    errors = np.loadtxt('./datafiles/errors%d.txt' %N)
    x = np.array(errors[:,0])
    abs_error = np.array(errors[:,1])
    rel_error = np.array(errors[:,2])
    x_log = np.log10(x)

    ax1.plot(x_log, abs_error, '-', color=colors[i-1], label='N=%d' %N)
    ax2.plot(x_log, rel_error, '-', color=colors[i-1], label='N=%d' %N)

ax1.set_title('Log of absolute errors for various number of steps N')
ax2.set_title('Log of relative errors for various number of steps N')

ax1.set(ylabel = 'Absolute error, log10(abs(u-v))')
ax2.set(ylabel = 'Relative error, log10(abs((u-v)/u))')

box = ax1.get_position()
ax1.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))

box = ax2.get_position()
ax2.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))

plt.xlabel('log(x)')
plt.savefig('./figures/error_plot_log.pdf')
plt.show()
plt.close()
