# Studying phase_transitions for various temperatures
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





# Phase transtions for L=40
dataL40 = np.loadtxt('./datafiles/phase_transitions_parallel_T(2.00-2.63)_L(40)_steps(20).txt', skiprows=1)
temperature_L40 = np.array(dataL40[:,0])
eps_L40 = np.array(dataL40[:,1])
m_L40 = np.array(dataL40[:,2])
cv_L40 = np.array(dataL40[:,3])
chi_L40 = np.array(dataL40[:,4])

# Phase transtions for L=60
dataL60 = np.loadtxt('./datafiles/phase_transitions_parallel_T(2.00-2.63)_L(60)_steps(20).txt', skiprows=1)
temperature_L60 = np.array(dataL60[:,0])
eps_L60 = np.array(dataL60[:,1])
m_L60 = np.array(dataL60[:,2])
cv_L60 = np.array(dataL60[:,3])
chi_L60 = np.array(dataL60[:,4])

# Phase transtions for L=80
dataL80 = np.loadtxt('./datafiles/phase_transitions_parallel_T(2.00-2.63)_L(80)_steps(20).txt', skiprows=1)
temperature_L80 = np.array(dataL80[:,0])
eps_L80 = np.array(dataL80[:,1])
m_L80 = np.array(dataL80[:,2])
cv_L80 = np.array(dataL80[:,3])
chi_L80 = np.array(dataL80[:,4])

# Phase transtions for L=100
dataL100 = np.loadtxt('./datafiles/phase_transitions_parallel_T(2.00-2.63)_L(100)_steps(20).txt', skiprows=1)
temperature_L100 = np.array(dataL100[:,0])
eps_L100 = np.array(dataL100[:,1])
m_L100 = np.array(dataL100[:,2])
cv_L100 = np.array(dataL100[:,3])
chi_L100 = np.array(dataL100[:,4])


# Plot epsilon
fig, axes = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.95,bottom=0.15,left=0.2,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$[T]=J/k_B$',ylabel='[$<\epsilon >$]=J')
plt.plot(temperature_L40, eps_L40, label='L=40',linewidth=1)
plt.plot(temperature_L60, eps_L60, label='L=60',linewidth=1)
plt.plot(temperature_L80, eps_L80, label='L=80',linewidth=1)
plt.plot(temperature_L100, eps_L100, label='L=100',linewidth=1)
plt.legend()
plt.grid()
plt.savefig("./figures/phase_trantions_eps.pdf")

# Plot magnetisation
fig, axes = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.95,bottom=0.15,left=0.15,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$[T]=J/k_B$',ylabel='[$<|m|>$]=1')
plt.plot(temperature_L40, m_L40, label='L=40',linewidth=1)
plt.plot(temperature_L60, m_L60, label='L=60',linewidth=1)
plt.plot(temperature_L80, m_L80, label='L=80',linewidth=1)
plt.plot(temperature_L100, m_L100, label='L=100',linewidth=1)
plt.legend()
plt.grid()
plt.savefig("./figures/phase_trantions_m.pdf")


# Plot heat capacity
fig, axes = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.95,bottom=0.15,left=0.15,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$[T]=J/k_B$',ylabel='[$C_v$]=$k_{B}^{2}$')
plt.plot(temperature_L40, cv_L40, label='L=40',linewidth=1)
plt.plot(temperature_L60, cv_L60, label='L=60',linewidth=1)
plt.plot(temperature_L80, cv_L80, label='L=80',linewidth=1)
plt.plot(temperature_L100, cv_L100, label='L=100',linewidth=1)
plt.legend()
plt.grid()
plt.savefig("./figures/phase_trantions_HC.pdf")

# Plot susceptibility
fig, axes = plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.95,bottom=0.15,left=0.15,right=0.95,hspace=0.2,wspace=0.2)
axes.set(xlabel='$[T]=J/k_B$',ylabel='[$\chi$]=$k_{B}/J$')
plt.plot(temperature_L40, chi_L40, label='L=40',linewidth=1)
plt.plot(temperature_L60, chi_L60, label='L=60',linewidth=1)
plt.plot(temperature_L80, chi_L80, label='L=80',linewidth=1)
plt.plot(temperature_L100, chi_L100, label='L=100',linewidth=1)
plt.legend()
plt.grid()
plt.savefig("./figures/phase_trantions_sus.pdf")

