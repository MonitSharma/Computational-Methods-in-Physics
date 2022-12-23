import matplotlib.pyplot as plt
import numpy as np


for i in range(1, 3):
    n = 10**i # n = 10, 100
    data = np.loadtxt('./datafiles/output%d.txt' %n)
    xhat = np.array(data[:,0])
    eigvec1 = np.array(data[:,1])
    eigvec2 = np.array(data[:,2])
    eigvec3 = np.array(data[:,3])
    analytical_eigvec1 = np.array(data[:,4])
    analytical_eigvec2 = np.array(data[:,5])
    analytical_eigvec3 = np.array(data[:,6])

    plt.plot(xhat, eigvec1, color='r', label='1 appr.')
    plt.plot(xhat, eigvec2, color='g', label='2 appr.')
    plt.plot(xhat, eigvec3, color='b', label='3 appr.')
    plt.plot(xhat, analytical_eigvec1, color='c', label='1 analyt.')
    plt.plot(xhat, analytical_eigvec2, color='m', label='2 analyt.')
    plt.plot(xhat, analytical_eigvec3, color='k', label='3 analyt.')

    plt.title('Eigenvectors of sym.tridiag matrix A with n=%d steps' %n)
    plt.ylabel('Eigenvectors, v(xhat)')
    plt.xlabel('xhat')
    plt.savefig('./figures/x_vs_eigvecn=%d.pdf' %n)
    plt.grid()
    plt.legend(bbox_to_anchor=(1.1, 1.05))
    plt.show()
