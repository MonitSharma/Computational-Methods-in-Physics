# relative error  agianst time for 5 different step sizes for 1 particle
# also plot the error convergence rate
import numpy as np
import matplotlib.pyplot as plt

SMALL_SIZE = 15
MEDIUM_SIZE = 15
BIGGER_SIZE = 15

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=13)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

max_error_EC = []
max_error_RK4 = []
stepsize =[]

for j in range(0,2):
    if (j == 0):
        method = 'EC'
    else:
        method = 'RK4'


    for i in range(1,6):
        iterations = 10**i



        data = np.loadtxt('./datafiles/%s_i_%d_d_100_p_1_pi_1_outputs_txyzv_pert_0_rs_0_f_0.0_w_v_0.0.txt'%(method,iterations), skiprows=1)

        t = np.array(data[:,0])
        x = np.array(data[:,1])
        y = np.array(data[:,2])
        z = np.array(data[:,3])
        v_x = np.array(data[:,4])
        v_y = np.array(data[:,5])
        v_z = np.array(data[:,6])

        if j==0:
            stepsize.append(t[-1]/iterations)
        #r = np.array(x,y,z)

        # --------- Error convergence rate (9.6)--------
        # max_error[i-1] = np.max(....)


        #------Analytical solution-----
        #Fetches the initial values and propertied for the particle
        x_0 = x[0]
        y_0 = 0
        z_0 = z[0]
        v_0 = v_y[0]

        m = 40.078
        q=1     #charge and mass

        #Defined constants from task

        B_0 = 9.65*10**1
        V_0=9.65 *10**8
        d=10**4

        #Constants
        omega_0 =(q*B_0)/m
        omega_z = np.sqrt((2*q*V_0)/(m*d**2))

        #constants
        omega_plus = (omega_0 + np.sqrt(omega_0**2 - 2*omega_z**2))/2
        omega_minus = (omega_0 - np.sqrt(omega_0**2 - 2*omega_z**2))/2

        #constants
        A_plus = (x_0*omega_minus + v_0)/(omega_minus- omega_plus)
        A_minus = -((v_0+x_0*omega_plus)/(omega_minus- omega_plus))

        #calculates exact values based on t(at timestep) and initial conditions
        x_exact = A_plus*np.cos(omega_plus*t) + A_minus*np.cos(omega_minus*t)
        y_exact = -(A_plus*np.sin(omega_plus*t) + A_minus*np.sin(omega_minus*t))
        z_exact = z_0*np.cos(omega_z*t)


        relative_error = np.sqrt((x-x_exact)**2+(y-y_exact)**2+(z-z_exact)**2)/np.sqrt((x_exact)**2+(y_exact)**2+(z_exact)**2)


        #appending delta max
        if method == 'EC':
            max_error_EC.append(max(np.sqrt((x-x_exact)**2+(y-y_exact)**2+(z-z_exact)**2)))

        if method == 'RK4':
            max_error_RK4.append(max(np.sqrt((x-x_exact)**2+(y-y_exact)**2+(z-z_exact)**2)))

        this_stepsize =t[-1]/iterations
        plt.yscale("log")
        plt.plot(t,relative_error, label ="Stepsize $h_k$= "+str(this_stepsize)+"$\mu s$")

    plt.xlabel("Time in $\mu s$")
    plt.ylabel("Relative error")
    plt.title(method)
    plt.legend(loc='lower right', ncol=2, fancybox=True, shadow=True, prop={'size':11})
    plt.grid()


    plt.subplots_adjust(
    top=0.93,
    bottom=0.135,
    left=0.16,
    right=0.985,
    hspace=0.2,
    wspace=0.2
    )
    plt.savefig('./figures/relative_error_%s.pdf'%(method), bbox_inches='tight')
    plt.clf()
    #plt.show()



r_err_sum_RK4 = 0
r_err_sum_EC = 0
for i in range(1,5):
    r_err_sum_RK4 += np.log(max_error_RK4[i]/max_error_RK4[i-1])/np.log(stepsize[i]/stepsize[i-1])
    r_err_sum_EC += np.log(max_error_EC[i]/max_error_EC[i-1])/np.log(stepsize[i]/stepsize[i-1])

r_err_sum_EC = r_err_sum_EC*(1/4)
r_err_sum_RK4 = r_err_sum_RK4 *(1/4)
print("Convergence rate for Euler is = ", r_err_sum_EC)
print("Convergence rate for RK4 is = ", r_err_sum_RK4)
