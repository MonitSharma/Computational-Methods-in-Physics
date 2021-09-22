from math import *
import matplotlib.pyplot as plt
array1 = []
array2 = []
array3 = []
def milne_simpson(f, f_actual, t, y, p, h):
	n = int((p-t)/h)
	time = [t+i*h for i in range(n+1)]

	y_values = [0]*(n+1)
	y_values[0] = y

	for i in range(1, 4):
		y_values[i] = rk4(f, time[i-1], y_values[i-1])

	for i in range(4, n+1):
		milne = y_values[i-4] + 4/3*h*(2*f(time[i-1], y_values[i-1]) - f(time[i-2], y_values[i-2]) + 2*f(time[i-3], y_values[i-3]))
		y_values[i] = y_values[i-2] + h/3*(f(time[i],milne) + 4*f(time[i-1], y_values[i-1]) + f(time[i-2], y_values[i-2]))

	print(f"Milne-Simpson Method")
	print("-"*70)
	print(f't\t\t|\tEstimate\t\t|\tActual')
	print("-"*70)
    
	for i in range(len(y_values)):
		print(f'{time[i]:.3f}\t\t|\t{y_values[i]:.9f}\t\t|\t{f_actual(time[i]):.9f}')
        array1.append(time)
        array2.append(y_values)
        array3.append(f_actual(time))
	print("-"*70)

def rk4(f, t, y):
	k1 = h*f(t, y)
	k2 = h*f(t + h/2, y + k1/2)
	k3 = h*f(t + h/2, y + k2/2)
	k4 = h*f(t + h, y + k3)
	return y + (k1 + 2*k2 + 2*k3 + k4)/6

if __name__ == '__main__':

	f = lambda t, y: (t**(-2))*(t*y - (y**2))

	f_actual = lambda t: t* (1/2 + log(t))**(-1)

	t = float(input("t0: "))
	y = float(input("y0: "))
	p = float(input("Evaluation point: "))
	h = float(input("Step size: "))

	milne_simpson(f, f_actual, t, y, p, h)


plt.plot(array1,array2)
plt.show()