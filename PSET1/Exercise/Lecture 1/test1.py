import matplotlib.pyplot as plt     

x = [1.000000,1.153049,1.251446,1.369597,1.517908,1.707702,1.956983,2.293680,2.763198,3.000000]
y = [2.000000,1.899793,1.816786,1.758056,1.721263,1.708262,1.722806,1.771618,1.865876,1.919219]



plt.title("Runge-Kutta Function Analysis")          # set the title of the graph
plt.xlabel("x")                                             # set the x label on the graph
plt.ylabel("y")                                             # set the y label on the graph
plt.plot(x, y, 'r--', label = "Runge Kutta")             # set the runge-kutta to be blue and label it
plt.legend()                                                # shows the legend on the graph
#plt.show()   


plt.savefig('Runge_kutta.png')