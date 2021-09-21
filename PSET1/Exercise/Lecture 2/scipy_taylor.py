from scipy.integrate import odeint
import numpy as np
from matplotlib import pyplot as plt

def f(t,x):
    
    return np.exp(np.power(-t,2))


x0 = 0
t = np.arange(0,2,0.01)

xs = odeint(f,x0,t)

plt.plot(t,xs,'-')
plt.plot(t,xs,'ro')
plt.xlabel('x values')
plt.ylabel('y values')
plt.title('Differential Equation')
plt.savefig('scipy_int.png')