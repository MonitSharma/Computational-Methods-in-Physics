'''
We are trying to solve an Ordinary Differential Equation
using the Euler's method. For this example we are taking
the example of:

dy/dt = -ky : IVP => y = y0 for t = 0

We will solve this without using the scipy
odeint package
'''

# import the basic libraries
import math
import numpy as np                # used for numerical analysis
import matplotlib.pyplot as plt   # used for plotting Graphs

# initialize the problem
x0 = -np.pi/2
t =0
x = x0
#k ,t , y =  1, 0.0, 5.0
#y0 = y
tt, yy = [], []
tf = 2                         # final value for the time

dt = 0.1                         # step size
#exact_sol = []
# describing the equation in function

def f(t,x):
    return -math.sin(t+x)


# now thw Euler method
exact_sol = []
while t<= tf:
    tt.append(t)
    yy.append(x)
    x = x+dt*f(t,x)
    exact = -math.asin((1-t**2)/(1+t**2)) - t
    exact_sol.append(exact)
    #exact = -math.asin((1-t**2)/(1+t**2))
    #exact_sol.append(exact)

    t += dt


# finding the Exact Result


#time 
exact = -math.asin((1-t**2)/(1+t**2)) - t

    
#print(exact)
# now the Plotting part
print(len(tt))

plt.plot(tt,yy,'o', label="dt=%.4f"%(dt))
plt.plot(tt, exact_sol,'k', label="Exact Solution")
plt.xlabel("time")
plt.ylabel('y')
plt.legend(loc='best')


# plotting the absolute error

plt.legend(loc='best')

plt.title("ODE by Euler Method")
plt.savefig('decay-euler1.png')