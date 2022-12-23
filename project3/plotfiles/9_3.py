# Phase space plots - position vs velocity in each direction for two particles!!!
import numpy as np
import matplotlib.pyplot as plt
import sys

# $axis $interaction $particle

duration = int(sys.argv[1])

# WITHOUT INTERACTIONS
data0x = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_0_outputs_xv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)
data0y = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_0_outputs_yv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)
data0z = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_0_outputs_zv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)

# ------Particle 1--------
# positions
x01 = np.array(data0x[:,0])
y01 = np.array(data0y[:,0])
z01 = np.array(data0z[:,0])
# velocities
vx01 = np.array(data0x[:,1])
vy01 = np.array(data0y[:,1])
vz01 = np.array(data0z[:,1])

# ------- Particle 2 ---------
# positions
x02 = np.array(data0x[:,2])
y02 = np.array(data0y[:,2])
z02 = np.array(data0z[:,2])
# velocities
vx02 = np.array(data0x[:,3])
vy02 = np.array(data0y[:,3])
vz02 = np.array(data0z[:,3])

# WITH INTERACTIONS

data1x = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_1_outputs_xv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)
data1y = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_1_outputs_yv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)
data1z = np.loadtxt('./datafiles/RK4_i_10000_d_%d_p_2_pi_1_outputs_zv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%duration, skiprows = 1)

# -------Particle 1 --------
# Positions
x11 = np.array(data1x[:,0])
y11 = np.array(data1y[:,0])
z11 = np.array(data1z[:,0])
# Velocities
vx11 = np.array(data1x[:,1])
vy11 = np.array(data1y[:,1])
vz11 = np.array(data1z[:,1])

# -------Particle 2 --------
# Positions
x12 = np.array(data1x[:,2])
y12 = np.array(data1y[:,2])
z12 = np.array(data1z[:,2])
# Velocities
vx12 = np.array(data1x[:,3])
vy12 = np.array(data1y[:,3])
vz12 = np.array(data1z[:,3])

# --------- PLOT w/subplots

# activate pylustrator
#pylustrator.start()

SMALL_SIZE = 9
MEDIUM_SIZE = 13
BIGGER_SIZE = 13

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


fig, axes = plt.subplots(2,3)#, sharey=True, sharex = True)
axes[0,0].set(ylabel='$Velocity(w/o interaction)$')#, title='No interaction')
axes[0,0].plot(x01, vx01, label='p1',linewidth=0.5)
axes[0,0].plot(x02, vx02, label='p2',linewidth=0.5)
#axes[0,0].ticklabel_format(axis="both", style="sci", scilimits=(0,0))
#axes[0,0].set_xlim([-3*10**3, 3*10**3])

axes[1,0].set(xlabel = 'x(t), [x(t)] = $ \mu m$', ylabel='$Velocity(w/interaction)$')#, title='Interaction')
axes[1,0].plot(x11, vx11, label='p1',linewidth=0.5)
axes[1,0].plot(x12, vx12, label='p2',linewidth=0.5)
axes[1,0].ticklabel_format(axis="both", style="sci", scilimits=(0,0))
#axes[1,0].set_xlim([-10**4, 10**4])

# y w/o
axes[0,1].set()#, title='No interaction y')
axes[0,1].plot(y01, vy01, label='p1',linewidth=0.5)
axes[0,1].plot(y02, vy02, label='p2',linewidth=0.5)
axes[0,1].ticklabel_format(axis="both", style="sci", scilimits=(0,0))
#axes[0,1].set_xlim([-3*10**3, 3*10**3])

# y w/
axes[1,1].set(xlabel = 'y(t), [y(t)] = $ \mu m$')
axes[1,1].plot(y11, vy11, label='p1 w/',linewidth=0.5)
axes[1,1].plot(y12, vy12, label='p2 w/',linewidth=0.5)
#axes[1,1].set_xlim([-10**4, 10**4])

# z w/o
axes[0,2].set()#, title='No interaction y')
axes[0,2].plot(z01, vz01, label='p1 w/o ',linewidth=0.5)
axes[0,2].plot(z02, vz02, label='p2 w/o',linewidth=0.5)
#axes[0,2].set_xlim([-3*10**3, 3*10**3])

# z w/
axes[1,2].set(xlabel = 'z(t), [z(t)] = $ \mu m$')#, title='interaction y')
axes[1,2].plot(z11, vz11, label='p1 w/',linewidth=0.5)
axes[1,2].plot(z12, vz12, label='p2 w/',linewidth=0.5)
#axes[1,2].set_xlim([-10**4, 10**4])


for axis in axes:
    for ax in axis:
        ax.ticklabel_format(axis="both", style="sci", scilimits=(0,0))

        #ax.set_xlim([-3**4, 3**4])

plt.subplots_adjust(
    top=0.9,
    bottom=0.1,
    left=0.1,
    right=0.95,
    hspace=0.3,
    wspace=0.3
    )

plt.savefig('./figures/phase_space%d.pdf'%duration)
#plt.show()
plt.close()
