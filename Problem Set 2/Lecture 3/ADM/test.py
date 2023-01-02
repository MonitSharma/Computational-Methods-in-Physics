from numpy import *
import numpy as np
from math import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


'''
Write down the function that we need to solve for 
x' = t^{-2}(tx-x^2)
'''

def f(t,x):
    return (t-3.2)*x + 8*t*exp((t-3.2)**2/2)*cos(4*t**2)

'''
another derivative of the function function
''' 

def dfx(t,x):
    return t-3.2



# analytical solution

def sol(t,t0,x0):
    
    return exp((t-3.2)**2/2)*(sin(4*t**2))


# backward Euler 
def beuler(t0, tn, n, x0):
    h = abs(tn - t0)/n
    t = linspace(0, tn, n+1)
    x = zeros(n+1)
    x[0] = x0

    for k in range(0,n):
        err = 1
        zold = x[k] + h*f(t[k], x[k])
        I = 0

        while err > 10**(-10) and I < 5:
            F = x[k] + h*f(t[k+1], zold) - zold
            dF = h*dfx(t[k+1], zold)-1
            znew = zold - F/dF
            err = abs(znew- zold)
            I += 1
        x[k+1] = znew
    return x




def RK4(t0,tn,n,x0):
    h = abs(tn-t0)/n
    t = linspace(t0, tn, n+1)
    x = zeros(n+1)
    x[0] = x0
    for i in range(0,n):
        K1 = f(t[i], x[i])
        K2 = f(t[i]+h/2, x[i]+K1*h/2)
        K3 = f(t[i] + h/2, x[i]+K2*h/2)
        K4 = f(t[i] + h, x[i]+K3*h)

        x[i+1] = x[i] + h*(K1+2*K2+ 2*K3+K4)/6
    return x



def AdBash3(t0,tn,n,x0):
    h = abs(tn-t0)/n
    t = linspace(t0, tn, n+1)
    x = zeros(n+1)

    x[0:3] = RK4(t0, t0+2*h,2,x0)
    K1 = f(t[1], x[1])
    K2 = f(t[0], x[0])

    for i in range(2,n):
        K3= K2
        K2= K1
        K1 = f(t[i], x[i])

        # prediction_resul
        x[i+1] = x[i] + h*(23*K1 - 16*K2+5*K3)/12
    return x





def PreCorr3(t0,tn,n,x0):
    h = abs(tn-t0)/n
    t = linspace(t0, tn, n+1)
    x = zeros(n+1)

    x[0:3] = RK4(t0, t0+2*h,2,x0)
    K1 = f(t[1], x[1])
    K2 = f(t[0], x[0])

    for i in range(2,n):
        K3= K2
        K2= K1
        K1 = f(t[i], x[i])

        # prediction_resul
        x[i+1] = x[i] + h*(23*K1 - 16*K2+5*K3)/12
    

        K0 = f(t[i+1],x[i+1])

        # correctorical

        x[i+1] = x[i] + h*(9*K0 + 19*K1 - 5*K2 + K3)/24

    return x


# print

fg = 1
n= 300
t0 = 0
tn = 6 
x0 = 0.75    
t = linspace(t0, tn , n+1)
xe = AdBash3(t0,tn, n, x0)
xb = beuler(t0,tn,n,x0)
xpc = PreCorr3(t0, tn, n,x0)

#figure(fg)

plot(t,xe,'red',label ='Adam Bashforth 3')
plot(t, xb, 'black', label ='Backward Euler')
plot(t,xpc,'cyan',label='Predictor Corrector')

t = linspace(t0, tn, 401)
xsol = sol(t,t0,x0)

plot(t,xsol, color='green', label='Exact')

#title('n = %d' %n)
#axis([0,tn, -60,40])
legend(loc='lower left')
plt.savefig('Test4.png')