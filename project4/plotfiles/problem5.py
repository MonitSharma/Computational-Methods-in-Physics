# Studying the burn-in time for L = 20
import numpy as np
#from plotlib import makeplots
import matplotlib.pyplot as plt


#Adjusting text size
SMALL_SIZE = 13
MEDIUM_SIZE = 15
BIGGER_SIZE = 17

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title


#fig, axes = plt.subplots(figsize = (1, 2))



# T = 1.0 ordered
dataT1o = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_1.0_ordered.txt', skiprows=1)
cycles_T1o = np.array(dataT1o[:,0])
eps_T1o = np.array(dataT1o[:,3])
m_T1o = np.array(dataT1o[:,4])

# T = 1.0 unordered
dataT1u = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_1.0_unordered.txt', skiprows=1)
cycles_T1u = np.array(dataT1u[:,0])
eps_T1u = np.array(dataT1u[:,3])
m_T1u = np.array(dataT1u[:,4])

# Plot T = 1.0 epsilon
fig, axes = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(top=0.90,bottom=0.15,left=0.2,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$n_{cycles}$',ylabel='$[<\epsilon >]$ = J')
plt.plot(cycles_T1o, eps_T1o, label='ordered',linewidth=1)
plt.plot(cycles_T1u, eps_T1u, label='unordered',linewidth=1)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()
plt.savefig('figures/burnin_eps_T1.pdf')

# Plot T = 1.0 magnetisation
fig, axes = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(top=0.90,bottom=0.15,left=0.2,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$n_{cycles}$',ylabel='$[<|m|>]$ = 1')
plt.plot(cycles_T1o, m_T1o, label='ordered',linewidth=1)
plt.plot(cycles_T1u, m_T1u, label='unordered',linewidth=1)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()
plt.savefig('figures/burnin_m_T1.pdf')

# T = 2.4 ordered
dataT2o = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_2.4_ordered.txt', skiprows=1)
cycles_T2o = np.array(dataT2o[:,0])
eps_T2o = np.array(dataT2o[:,3])
m_T2o = np.array(dataT2o[:,4])

# T = 2.4 unordered
dataT2u = np.loadtxt('./datafiles/ncyc_1e4_L_20_T_2.4_unordered.txt', skiprows=1)
cycles_T2u = np.array(dataT2u[:,0])
eps_T2u = np.array(dataT2u[:,3])
m_T2u = np.array(dataT2u[:,4])

#axes[0,1].set(ylabel='$<\epsilon >$') #, title='No interaction')
fig, axes = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(top=0.90,bottom=0.15,left=0.2,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$n_{cycles}$',ylabel='$[<\epsilon >$] = J')
plt.plot(cycles_T2o, eps_T2o, label='ordered',linewidth=1)
plt.plot(cycles_T2u, eps_T2u, label='unordered',linewidth=1)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()
plt.savefig('figures/burnin_eps_T24.pdf')

# Plot T = 1.0 magnetisation
#axes[1,0].set(ylabel='$<m >$') #, title='No interaction')
fig, axes = plt.subplots(figsize = (6, 5))
plt.subplots_adjust(top=0.90,bottom=0.15,left=0.2,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$n_{cycles}$',ylabel='$[<|m|>]$ =1')
plt.plot(cycles_T1o, m_T2o, label='ordered',linewidth=1)
plt.plot(cycles_T1u, m_T2u, label='unordered',linewidth=1)
plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
plt.grid()
plt.legend()
plt.savefig('figures/burnin_m_T24.pdf')
