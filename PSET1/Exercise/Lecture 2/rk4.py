import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# define equations

def dydt(t, y): return np.sin(t)**2*y
def f(t): return 2*np.exp(0.5*(t-np.sin(t)*np.cos(t)))


def RK3(t, y, h):
    # compute approximations
    k_1 = dydt(t, y)
    k_2 = dydt(t+h/2, y+(h/2)*k_1)
    k_3 = dydt(t+h/2, y+h*(-k_1 + 2*k_2))

    # calculate new y estimate
    y = y + h * (1/6) * (k_1 + 4 * k_2 + k_3)
    return y


def RK4(t, y, h):
    # compute approximations
    k_1 = dydt(t, y)
    k_2 = dydt(t+h/2, y+(h/2)*k_1)
    k_3 = dydt(t+h/2, y+(h/2)*k_2)
    k_4 = dydt(t+h, y+h*k_3)

    # calculate new y estimate
    y = y + h * (1/6)*(k_1 + 2*k_2 + 2*k_3 + k_4)
    return y


# Initialization
ti = 0
tf = 5
h = 0.5
n = int((tf-ti)/h)
t = 0
y = 2
y_rk3 = 2

print("t \t\t yRK3 \t\t yRK4 \t\t f(t)")
print(f"{t:.1f} \t\t {y_rk3:4f} \t\t {y:4f} \t\t {f(t):.4f}")

t_plot = []
y_RK3 = []
y_RK4 = []
y_analytical = []


##
for i in range(1, n+1):
    t_plot.append(t)
    y_RK4.append(y)
    y_RK3.append(y_rk3)
    y_analytical.append(f(t))

    # calculate new y estimate
    y = RK4(t, y, h)
    y_rk3 = RK3(t, y_rk3, h)

    t += h
    print(f"{t:.1f} \t\t {y_rk3:4f} \t\t {y:4f} \t\t {f(t):.4f}")


t_plot.append(t)
y_RK3.append(y_rk3)
y_RK4.append(y)
y_analytical.append(f(t))


# Visualization
fig, (ax, ax2) = plt.subplots(2, 1)
ax.plot(t_plot, y_analytical, 'ro', label='Analytical solution')
ax.plot(t_plot, y_RK4, '.b', label='Fourth-order Runge-Kutta estimate')
ax.plot(t_plot, y_RK3, '.g', label='Third-order Runge-Kutta estimate')
ax.set_ylabel("y", fontsize=18)
ax.legend()

ax2.plot(t_plot, np.abs(np.array(y_analytical) -
                        np.array(y_RK4)), '.-b', label='Fourth-order Runge-Kutta')
ax2.plot(t_plot, np.abs(np.array(y_analytical) -
                        np.array(y_RK3)), '.-g', label='Third-order Runge-Kutta')
ax2.set_ylabel("Abs Error", fontsize=18)
ax2.legend()
ax2.set_xlabel("t", fontsize=18)
plt.savefig('runge_kutta_analytical.png', dpi=300, bbox_inches='tight')