import numpy as np
import matplotlib.pyplot as plt
from scipy import stats




lattice_sizes = [40,60,80,100]

T_c =[]
L_inverse=[]

for L_size in lattice_sizes:

	name_template = "phase_transitions_parallel_T(2.00-2.63)_L("+str(L_size)+")_steps(20).txt"

	data = np.loadtxt('./datafiles/'+name_template, skiprows=1)
	T = np.array(data[:,0])
	chi = np.array(data[:,4])

	max_chi = np.amax(chi) 				#Finds max value of chi
	max_index = np.where(chi==max_chi)	#Finds what the index of the max value is

	T_c.append(T[max_index][0])			#Append only value
	L_inverse.append(1/L_size)

	#print('For L=',L_size,'found max sus=',max_chi,', at temperature T_c=',T[max_index][0])



slope, intercept, r, p, intercept_err = stats.linregress(L_inverse,T_c)

T_c.insert(0,intercept)
L_inverse.insert(0,0)

T_c = np.array(T_c)
L_inverse = np.array(L_inverse)



#print(result.intercept)
#print(T_c)
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


plt.subplots(figsize=(6,5))
plt.subplots_adjust(top=0.95,bottom=0.15,left=0.15,right=0.95,hspace=0.2,wspace=0.2)
plt.plot(L_inverse[1:],T_c[1:],'ro',label='data points')
plt.plot(L_inverse[0],T_c[0],'bx',label='Intercept='+str(round(intercept,3))+'$\pm$'+str(round(intercept_err,3)))
plt.plot(L_inverse, (intercept + slope*L_inverse), 'g', label='fitted line')
plt.xlabel("1/L")
plt.ylabel("[$T_c$]= J/$k_B$")
plt.legend()
plt.grid()
#print('The estimate critical temperature for an infinite lattice is: ',result.intercept)
plt.savefig("./figures/critical_temp_inf.pdf")