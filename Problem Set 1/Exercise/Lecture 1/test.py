'''
Use the Taylor Series method of order 4 to compute $\int_0^2 e^-s^2 ds
by solving the Initial Value problem on the interval t \in [0,2]
Given dx/dt = e^-t^2
and   x(0) = 0
Check from the table of error function x(2) = (pi^(1/2)/2). erf(2) = 0.8820813907
We need to find the derivative till the fourth order of the Given function, and add them to a Taylor
Series to form the solution
'''

# importing the basic libraries
from sympy import *                             # used for symbols in mathematics
import math                                     # used for math function
import numpy as np                              # used for numerical analysis and other
import matplotlib.pyplot as plt                 # used for plotting the graphs and saving them



# we will begin by defining the function 
# given to us , and then the four orders
# of derivatives

'''
The first order derivative , and this function
is given to us
'''
def f1(t,x):
    return exp((-t**2))


# now we differentiate the same function and find the 
# next orders of derivatives
def f2(t,x):
    return -2*t*f1(t,x)

def f3(t,x):
    return 4*np.power(t,2)*f1(t,x) - 2*f1(t,x)

def f4(t,x):
    return 12*t*f1(t,x) - 8*np.power(t,3)*f1(t,x)


# also the function that finds the Error

def error(dt,m):
    return ((dt**m)/factorial(m+1)) * (xx[len(xx)-1] - xx[len(xx)-2])



# and the exact solution to the problem(for reference)
def exact(t):
    return math.sqrt(np.pi)*math.erf(t)*0.5



# Now we write down the values known to us
m = 4                       # order parameter
dt = 0.2                    # step Size
t,x = 0 ,0                    # initial values of both
tf = 2                      # final point to reach
xe = 0
# initialised arrays for graph plotting
tt = []
xx = []
xxe = []

# now the loop

while t<=tf:
    tt.append(t)                     # append the empty list for x-axis for each time step
    xx.append(x)                     # append this empty list for the y-axis for each value of x
    xxe.append(xe)
    x = x+dt*f1(t,x) + np.power(dt,2)/2*f2(t,x) + np.power(dt,3)/6*f3(t,x) + np.power(dt,4)/24*f4(t,x)       # the taylor series
    #xe = exact(t)
    t+= dt                           # increase the time step by 0.1





# Plotting and Visualization
import matplotlib.pyplot as plt
import numpy as np
import scipy

from scipy.special import erf
t = np.linspace(0,2,10)
#plt.plot(t,erf(t),color="#A60628", linewidth = 2, label ='erf(x)')
#plt.savefig("test.png")

plt.plot(tt,xx,'o',linewidth=2,color='red', label ='Taylor Method Order = 4')
plt.plot(t,math.sqrt(np.pi)*erf(t)*0.5,'b-',linewidth=2,label= 'Exact Solution')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.legend()
plt.savefig('Test1.png')



# checking for the value at t = 2
# writing all the values 

print("The Value that x(t) takes for different t:\n", xx)

# printing the error

print("The Error Value is:\n ",error(dt,m))

#  difference between exact value and Taylor approximated

print("The Difference between the Exact solution at t = 2 and Taylor approximated is:\n", math.erf(2)-xx[len(xx)-1])


