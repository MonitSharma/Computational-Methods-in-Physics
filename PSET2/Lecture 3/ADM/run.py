from numpy import *
import numpy as np
from math import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def f(t,x):
    return (t**(-2))*(t*x - (x**2))



n= 300
t0 = 1
tn = 3 
x0 = 2 
h = abs(tn - t0)/n
t = t0
#t = linspace(0, tn, n+1)
x = zeros(n+1)

#print(f(2,1))

#while (t<=tn):
#    print(f(t,x))
#    t = t+h



#t = linspace(t0, tn, n+1)


t = linspace(t0, tn, n+1)
print(t[2]-t[1])




#print(exact(3))


