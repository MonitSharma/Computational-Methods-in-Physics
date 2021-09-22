'''
Using the Shooting method to solve
y" = -3y
with the boundary conditions y(0) = 7 and y(2*pi) = 0
and we know the exact solution
''' 


# importing stuff
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


x = np.linspace(0,2*np.pi, 100)
y_exact = 7 * np.cos(np.sqrt(3) * x) - 7 * np.cos(2 * np.pi * np.sqrt(3) ) / np.sin(2 * np.pi * np.sqrt(3) )* np.sin(np.sqrt(3)*x)


# code our second order ODE to two first order ODE 

def equations(x,y):
    yprime = np.zeros(2)

    yprime[0] = y[1]
    yprime[1] = -3* y[0]

    return yprime

'''
We will use an iteration scheme to adjust the 
values of the initial conditions such that we
match the boundary conditions. What type of 
iteration scheme to use is up to you. I will 
define high and low bounds and set our natural 
guess As the meaning of those two points. I will
then adjust the bounds depending on whether we 
overshoot or undershoot our boundary. The exact
way to adjust the bounds is problem dependent.
'''

tol = 1e-6
max_iters = 100 
low = -10  
high = 10
count = 0 

while count <= max_iters:
    count = count + 1
    xspan = (x[0], x[-1])
    
    #  Use the midpoint between high and low as our guess
    yprime0 = np.mean([low, high])
    
    #  Set the initial condition vector to be passed into the solver
    y0 = [7, yprime0 ]

    # Solve the system using our guess
    sol = solve_ivp(equations, xspan, y0, t_eval = x)

    #  For ease of use, extract the function values from the solution object.
    y_num = sol.y[0, :]

    #  Check to see if we within our desire tolerance
    if np.abs(y_num[-1]) <= tol:
        break
    
    #  Adjust our bounds if we are not within tolerance
    if y_num[-1] < 0:
        high = yprime0
    else:
        low = yprime0
        
    #print(count, y_num[-1])
    
#  Plot the solution and compare it to the analytical form defined above
plt.plot(x, y_exact, 'b-', label='Exact')
plt.plot(x, y_num, 'ro', label='Numeric')
plt.plot([0, 2*np.pi], [7,0], 'ro')
plt.grid(True)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
#plt.show()
plt.savefig('shooting_1.png')