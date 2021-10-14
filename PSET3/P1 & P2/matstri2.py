# EqStringAnimateMat.py:  Animated leapfrog solution Vibrating string using MatPlotLib

from numpy import *
import numpy as np, matplotlib.pyplot as plt, matplotlib.animation as animation

rho = 0.01;   ten = 40.; c = sqrt(ten/rho)               # density, tension 
c1 = c;       ratio =  c*c/(c1*c1)                       # CFL criterium = 1
xi = np.zeros((101,3), float)                            # Declaration 
k = range(0,101)

def Initialize():                                       # Initial conditions
   for i in range(0, 81):    xi[i, 0] = 0.00125*i         
   for i in range (81, 101): xi[i, 0] = 0.1 - 0.005*(i - 80)                  # second part of string
 
def animate(num):                                
    for i in range(1, 100):              
       xi[i,2] = 2.*xi[i,1]-xi[i,0]+ratio*(xi[i+1,1]+xi[i-1,1]-2*xi[i,1])
    line.set_data(k,xi[k,2])                                # Data to plot ,x,y           
    for m in range (0,101):                                
       xi[m, 0] = xi[m, 1]                                    # Recycle array 
       xi[m, 1] = xi[m, 2]
    return line        
         
Initialize()                                             # Plot initial string   
fig = plt.figure()                            
ax = fig.add_subplot(111, autoscale_on=False, xlim=(0, 101), ylim=(-0.15, 0.15))
ax.grid()                                                        # Plot  grid
plt.title("Vibrating String")
line, = ax.plot(k, xi[k,0], lw=2)                 
for i in range(1,100): 
    xi[i,1] = xi[i,0] + 0.5*ratio*(xi[i+1,0] + xi[i-1,0] -2*xi[i,0])  
ani = animation.FuncAnimation(fig, animate,1)             # Dummy argument: 1       
plt.show()             
print("finished")