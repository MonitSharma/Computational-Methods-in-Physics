import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import matplotlib.animation as an

import networkx as nx

#Butcher Tableau
alpha = np.array([0.0, 1.0/4.0, 3.0/8.0, 12.0/13.0, 1.0, 1.0/2.0])
beta = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
		[1.0/4.0, 0.0, 0.0, 0.0, 0.0, 0.0],
		[3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0, 0.0],
		[1932.0/2197.0, (-7200.0)/2197.0, 7296.0/2197.0, 0.0, 0.0, 0.0],
		[439.0/216.0, -8.0, 3680.0/513.0, (-845.0)/4104.0, 0.0, 0.0],
		[(-8.0)/27.0, 2.0, (-3544.0)/2565.0, 1859.0/4104.0, (-11.0)/40.0, 0.0]])
c = np.array([25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, (-1.0)/5.0, 0.0]) # coefficients for 4th order method
c_star = np.array([16.0/135.0, 0.0, 6656.0/12825.0, 28561.0/56430.0, (-9.0)/50.0, 2.0/55.0]) # coefficients for 5th order method
cerr=c-c_star

def rkf(f,x,h,alpha,beta,cerr,c_star,*args,**kwargs):
    '''
    To perform integration stepwise using the Runge-Kutta-Fehlberg 4(5) method with adaptive step size control
    Parameters
    ----------
    f : function
        calculating the rhs of the ode
    x : numpy array
        initial state of the system
    h : float
        initial stepsize
    alpha : numpy array
        alpha values of the RK-4(5) method
    beta : numpy array
        beta values of the RK-4(5) method
    cerr : numpy array
        difference in the c cooefficients of the RK-4 and RK-5 methods
    c_star : numpy array
        c values of the RK-5 method
    *args : additional arguments to passdown to function f
    **kwargs : additional keywork arguments to pass down to function f
    Returns
    -------
    X : numpy array
        calculated values of the system for each time step
    T : numpy array
        time values corresponding to each time step in the system
    H : numpy array
        adaptive step size values for every calculated step
    EPS : numpy array
        error epsilon at each step(between the calculated values from RK-4 and RK-5 method)
    '''
    T = np.zeros(N_t)
    X = np.zeros((N,N_t))
    H = np.zeros(N_t)
    EPS = np.zeros(N_t)
    X[:,0] = x
    H[0] = h
    
    i = 0
    while i < N_t-1 : 
        k1 = f(X[:,i],*args,**kwargs)
        k2 = rhs(X[:,i]+h*beta[1,0]*k1, *args,**kwargs)
        k3 = rhs(X[:,i]+h*(beta[2,0]*k1 + beta[2,1]*k2),*args,**kwargs)
        k4 = rhs(X[:,i]+h*(beta[3,0]*k1 + beta[3,1]*k2 + beta[3,2]*k3),*args,**kwargs)
        k5 = rhs(X[:,i]+h*(beta[4,0]*k1 + beta[4,1]*k2 + beta[4,2]*k3 + beta[4,3]*k4),*args,**kwargs)
        k6 = rhs(X[:,i]+h*(beta[5,0]*k1 + beta[5,1]*k2 + beta[5,2]*k3 + beta[5,3]*k4 + beta[5,4]*k5),*args,**kwargs)
        errorfield = h*(cerr[0]*k1 + cerr[1]*k2 + cerr[2]*k3 + cerr[3]*k4 + cerr[4]*k5 + cerr[5]*k6)
        max_error = np.absolute(errorfield).max()
        
        if (max_error <= tol):
            X[:,i+1] = X[:,i] + h*(c_star[0]*k1 + c_star[1]*k2 + c_star[2]*k3 + c_star[3]*k4 + c_star[4]*k5 + c_star[5]*k6)
            h = safe * h* (tol/max_error)**0.2
            H[i+1]=h
            T[i+1]=T[i]+h
            EPS[i+1]=max_error
            i+=1
        else:
            h=safe*h*(tol/max_error)**0.25
    
    return X,T,H,EPS

def rhs(x,A,d,u,al,gm,b):
    '''
    TO calculate the RHS of the given ode
    Parameters
    ----------
    x : numpy array
        curent state of the system
    A : numpy array
        the adjacency matrix of the system graph
    d : float
        parameter, resistance to becoming opinionated
    u : float
        control parameter, social influence
    al : float
        parameter, self reinforcement
    gm : float
        parameter, cooperative/competitive
    b : float
        parameter, input bias
    Returns
    -------
    x_dot : numpy array
        RHS of the ode. The derivative of the current state.
    '''
    x1=np.zeros(np.shape(x))
    for i in range(np.shape(A)[0]):
        for j in range(np.shape(A)[1]):
            if (i!=j):
                x1[i]+=gm*A[i,j]*x[j]
            else:
                x1[j]+=al*x[j]
    x_dot=-d*x+u*np.tanh(x1)+b
    return x_dot


N_t = 200  #no. of steps taken in the adaptive step size
N = 10   #no. of nodes
tol = 0.0000001     #desired accuracy/tolerance
safe = 0.84     #safety factor

#Adjacency Matrix : 

A=np.zeros((N,N))

#TOPOLOGIES

# #Circle 
# for i in range(N-1):
#     A[i,i+1]=1
#     A[i+1,i]=1
    
#     A[N-1,0]=1
#     A[0,N-1]=1

# #Wheel
# for i in range(N-2):
#     A[i,i+1]=1
#     A[i+1,i]=1
#     A[i,N-1]=1
#     A[N-1,i]=1
# A[N-2,N-1]=1
# A[N-1,N-2]=1
# A[N-2,0]=1
# A[0,N-2]=1

# #Path
# for i in range(N-1):
#     A[i,i+1]=1
#     A[i+1,i]=1

#Star
for i in range(N-1):
    A[i,N-1]=1
    A[N-1,i]=1

#Plotting : 

h_in = 0.05 #initial stepsize
x0=(np.random.uniform(-1.0,1.0,size=N)) #inital random state of the system


x,t,h_val,err_val = rkf(rhs,x0,h_in,alpha,beta,cerr,c_star,A,d=0.5,u=0.26,al=1.2,gm=1.3,b=0.0)  
#Disagreement when gm(gamma) is negative and Agreement when gm(gamma) is positive

# for i in range(N):
#     plt.plot(t,x[i,:])
    
# plt.xlabel("time t")
# plt.ylabel("opinions x")
# plt.title("Star Topology Disagreement")
# #plt.savefig("Star_Topology_Disagreement.jpg")
# plt.show()

# plt.plot(t,h_val)
# plt.xlabel("time t")
# plt.ylabel("stepsize h")
# plt.title("Adaptive Step Size")
# plt.show()

# plt.plot(t,err_val)
# plt.xlabel("time t")
# plt.ylabel("error $\epsilon$")
# plt.title("Change in error")
# plt.show()

#Graph : 

rows,cols=np.where(A==1.)
edges=zip(rows.tolist(),cols.tolist())
G=nx.Graph()
G.add_edges_from(edges)
plot_positions=nx.drawing.spring_layout(G)

vmin=-1
vmax=1
norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
cmap=plt.get_cmap('coolwarm')
sm=plt.cm.ScalarMappable(cmap=cmap,norm=norm)
sm.set_array([])

#Animation. Save Video

fig=plt.gcf()

def animate(i):
    '''
    Function to iterate over in order to animate the graphs
    Parameters
    ----------
    i : int
        interative variable for the FuncAnimation() function to animate the graph
    Returns
    -------
    None.
    '''
    nx.draw(G,pos=plot_positions,node_size=500,node_color=x[:,i],cmap='coolwarm',vmin=vmin,vmax=vmax)


plt.colorbar(sm)
anim = an.FuncAnimation(fig, animate, frames=200, blit=False)
writervideo = an.FFMpegWriter(fps=10) 
anim.save('Star_Topology_Agreement.mp4', writer=writervideo)
