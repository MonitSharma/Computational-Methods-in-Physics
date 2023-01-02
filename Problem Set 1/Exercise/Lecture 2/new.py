import math
import matplotlib.pyplot as plt

t = 0
h = 0.1
tf = 2
exact_sol = []
while (t<=tf):
    exact = -math.asin((1-t**2)/(1+t**2)) - t
    exact_sol.append(exact)
    t = t+h


print(len(exact_sol))