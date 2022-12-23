# Plot of 2 particles and their motion in xy-plane with and without particle interactions
import numpy as np
import matplotlib.pyplot as plt
import sys

# pi:   0: no interactions, 1:interactions

duration = int(sys.argv[1])

data0 = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_0_outputs_xy_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows=1)

data1 = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_1_outputs_xy_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows=1)


# No interactions: $axis $interaction $particle
x01 = np.array(data0[:,0])
y01 = np.array(data0[:,1])
x02 = np.array(data0[:,2])
y02 = np.array(data0[:,3])

# Interactions: $axis $interaction $particle
x11 = np.array(data1[:,0])
y11 = np.array(data1[:,1])
x12 = np.array(data1[:,2])
y12 = np.array(data1[:,3])


SMALL_SIZE = 13
MEDIUM_SIZE = 17
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


fig, axes = plt.subplots(1,2, sharey=True,figsize=(8,4))#,figsize=(7, 7)) # default:6.4, 4.8 sharex = True, sharey=True

# Plotting no interaction vs interaction:

axes[0].plot(x01,y01,label='Particle 1',linewidth=0.5)
axes[0].plot(x02,y02,label='Particle 2',linewidth=0.5)

axes[1].plot(x11,y11,label='Particle 1',linewidth=0.5)
axes[1].plot(x12,y12,label='Particle 2',linewidth=0.5)


for ax in axes:
    ax.set_xlim([-0.22*10**4, 0.22*10**4])
    ax.set_ylim([-0.22*10**4, 0.22*10**4])
    #ax.add_patch(plt.Circle((0, 0), 10**4, linestyle="--", color='grey', fill=False)) # for penningtrap circle
    ax.set(adjustable='box', aspect='equal')
    ax.grid()
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

axes[0].set(xlabel = 'x(t), [x(t)] = $\mu m$', ylabel='y(t), [y(t)] = $\mu m$', title='No interaction')
axes[1].set(xlabel = 'x(t), [x(t)] = $\mu m$', title = 'Interaction')
#plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)#plt.legend()

#plt.legend(bbox_to_anchor=(1.05, 1), loc='best') # x,y position

#plt.legend(loc='best')

# handles, labels = ax.get_legend_handles_labels()
# fig.legend(handles, labels, loc='upper center')


plt.subplots_adjust(
    top=0.915,
    bottom=0.165,
    left=0.13,
    right=0.95,
    hspace=0.2,
    wspace=0.2
)

plt.savefig('./figures/xy_%d.pdf'%duration, bbox_inches='tight')
#plt.show()
