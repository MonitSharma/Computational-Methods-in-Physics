from sympy import *
from functools import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import subprocess
import copy

import sys



"""
Implementation of Numerical Methods
Center for Informatics, Federal University of Pernambuco (CIn/UFPE)
@author: Ermano A. Arruda <eaa3@cin.ufpe.br>
Implemented methods:
----One-step/Stepwise/Starting Methods----
1) Euler (error ~ h^2)
2) Improved Euler Method (Modified Euler) (error ~ h^3)
3) Backward Euler Method (error ~ h^2)
4) Runge-Kutta Method (error ~ h^4)
5) Three Term Taylor Series Method (error ~ h^3)
----Multistep or Continuing Methods----
6) Adams-Bashforth[1,2,3,4] (error = h^2, h^3, h^4, h^5 - respectively ) [TODO: Check error correctness]
7) Adams-Multon[1,2,3,4] (error = h^2, h^3, h^4, h^5 - respectively ) [TODO: Check error correctness]
8) Preditor-Corrector[1,2,3,4] (error = h^2, h^3, h^4, h^5 - respectively ) [TODO: Check error correctness]
9) Backward Differentiation (BackDiff[1,2,3,4]) (error = h^2, h^3, h^4, h^5 - respectively ) [TODO: Check error correctness]
"""

class DSolver:

    def __init__(self, yd_expression_str, phi_expr_str = None, episolon = 0.1):

        self.x_symb, self.y_symb, self.yd_symb = symbols("x y yd")
        self.yd_expression_str = yd_expression_str
        self.yd_expr = sympify(yd_expression_str)

        self.yd_func = lambdify((self.x_symb,self.y_symb),self.yd_expr,"numpy")

        self.x = np.zeros(1)
        self.y = np.zeros
        self.yd = np.zeros
        self.phi = np.zeros
        self.error = np.zeros
        self.accerror = np.zeros

        self.iyd = 0

        self.episolon = episolon

        self.phi_func = None
        if( phi_expr_str != None ):
            self.phi_expr = sympify(phi_expr_str)
            self.phi_func = lambdify(self.x_symb,self.phi_expr)

    # Initialization
    def __initialize__(self,x0,y0,n):
        
        self.x = self.y = self.yd = self.phi = self.error = self.accerror = None

        self.x = np.zeros(n+1)
        self.y = np.zeros(n+1)
        self.yd = np.zeros(n+1)
        self.phi = np.zeros(n+1)
        self.error = np.zeros(n+1)
        self.accerror = np.zeros(n+2)

        self.iyd = 0

        self.x[0] = x0
        self.y[0] = y0
        self.error[0] = 0

        if self.phi_func != None:
            self.phi[0] = self.phi_func(0)
            self.error[0] = self.phi[0] - self.y[0]
            self.accerror[0] = self.accerror[0]

    # Yd (derivative of y) -- Unused
    def __yd__(self,h,i):
        if self.iyd < i:
            self.iyd = i
            self.yd[i] = self.yd_func(i*h,self.y[i])*h
        
        return self.yd[i]


    # Default Euler method

    def __euler__(self,h,i):

        self.y[i] = self.y[i-1] + self.yd_func((i-1)*h,self.y[i-1])*h

        return self.y[i]

    # Modified (Improved) Euler method

    def __improved_euler__(self,h,i):
        self.y[i] = self.y[i-1] + h*(self.yd_func((i-1)*h,self.y[i-1])+self.yd_func(i*h,self.y[i-1] + self.y[i-1]*h))*0.5
        return self.y[i]


    # Backward Euler method

    def __backward_euler__(self,h,i):
        
        yi = symbols("yi")

        exp = self.y[i-1] + self.yd_expr.subs({self.y_symb:yi, self.x_symb:i*h})*h
        equ = Eq(exp,yi)
        
        self.y[i] = solve(equ,yi)[0]

        return self.y[i]
    

    # Runge-Kutta method

    def __runge_kutta__(self,h,i):
        kn1 = kn2 = kn3 = kn4 = 0
        
        kn1 = self.yd_func((i-1)*h,self.y[i-1])
        kn2 = self.yd_func(((i-1)+0.5)*h,self.y[i-1] + 0.5*h*kn1)
        kn3 = self.yd_func(((i-1)+0.5)*h,self.y[i-1] + 0.5*h*kn2)
        kn4 = self.yd_func(i*h,self.y[i-1] + h*kn3)

        self.y[i] = self.y[i-1] + h*(kn1 + 2*kn2 + 2*kn3 + kn4)/6

        return self.y[i]

    
    # Three term Taylor Series method

    def __3term_taylor_series__(self,f_x,f_y,h,i):


        yd = self.yd_func((i-1)*h,self.y[i-1])
        
        ydd = f_x((i-1)*h, yd) + f_y((i-1)*h, yd)*yd

        self.y[i] = self.y[i-1] + h*yd + ((h**2)/2)*ydd

        return self.y[i]


    # Three term Adams-Bashforth Series method

    def __adams_bashforth__(self,p, h,i):

        integral = 0

        if i <= p:
            self.__runge_kutta__(h,i)
        else:

            if p == 1:

                    ydn = self.yd_func((i-1)*h,self.y[i-1])
                    ydn_1 = self.yd_func((i-2)*h,self.y[i-2])

                    integral = (3*ydn*0.5 - ydn_1*0.5)*h

            elif p == 2:
                    ydn = self.yd_func((i-1)*h,self.y[i-1])
                    ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                    ydn_2 = self.yd_func((i-3)*h,self.y[i-3])

                    integral = h*((23.0/12.0)*ydn - (4.0/3.0)*ydn_1 + (5.0/12.0)*ydn_2)

            elif p == 3:

                    ydn = self.yd_func((i-1)*h,self.y[i-1])
                    ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                    ydn_2 = self.yd_func((i-3)*h,self.y[i-3])
                    ydn_3 = self.yd_func((i-4)*h,self.y[i-4])

                    integral = h*((55*ydn - 59*ydn_1 + 37*ydn_2 - 9*ydn_3)/24.0)

            elif p == 4:
                    ydn = self.yd_func((i-1)*h,self.y[i-1])
                    ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                    ydn_2 = self.yd_func((i-3)*h,self.y[i-3])
                    ydn_3 = self.yd_func((i-4)*h,self.y[i-4])
                    ydn_4 = self.yd_func((i-5)*h,self.y[i-5])

                    integral = h*((1901/720)*ydn - (1387/360)*ydn_1 + (109/30)*ydn_2 - (637/360)*ydn_3 + (251/720)*ydn_4)


            self.y[i] = self.y[i-1] + integral

        return self.y[i]

    # Three term Adams-Multon Series method

    def __adams_multon__(self,p, h,i):

        integral = 0


        if i < p:
            self.__runge_kutta__(h,i)
            return self.y[i]



        yi = symbols("yi")

        ydn1 = self.yd_expr.subs({self.y_symb:yi, self.x_symb:i*h})


        #if polinomial degree 1
        if p == 1:

                ydn = self.yd_func((i-1)*h,self.y[i-1])

                integral = (ydn1 + ydn)*0.5*h

        #if polinomial degree 2
        elif p == 2:

                ydn = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])

                integral = h*((5.0*ydn1 + 8.0*ydn - ydn_1)/12.0)
        #if polinomial degree 3
        elif p == 3:

                ydn = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                ydn_2 = self.yd_func((i-3)*h,self.y[i-3])

                integral = h*((9.0*ydn1 + 19.0*ydn - 5.0*ydn_1 + ydn_2)/24.0)
        #if polinomial degree 4
        elif p == 4:

                ydn  = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                ydn_2 = self.yd_func((i-3)*h,self.y[i-3])
                ydn_3 = self.yd_func((i-4)*h,self.y[i-4])

                integral = ((251.0*ydn1)/720.0 + (646.0*ydn)/720.0 - (264.0*ydn_1)/720.0 + (106.0*ydn_2)/720.0 - (19.0*ydn_3)/720.0)*h


        exp = self.y[i-1] + integral
        equ = Eq(exp,yi)

        self.y[i] = solve(equ,yi)[0]

        return self.y[i]


    # Prediction correction
    def __predictor_corrector__(self,p, h,i):

        integral = 0

        # Prediction step
        self.__adams_bashforth__(p,h,i)

        


        if i < p:
            return self.y[i]


        #Correction step

        #if polinomial degree 1
        if p == 1:

                ydn1 = self.yd_func(i*h,self.y[i])
                ydn = self.yd_func((i-1)*h,self.y[i-1])

                integral = (ydn1 + ydn)*0.5*h

        #if polinomial degree 2
        elif p == 2:

                ydn1 = self.yd_func(i*h,self.y[i])
                ydn = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])

                integral = h*((5*ydn1 + 8*ydn - ydn_1)/12.0)
        #if polinomial degree 3
        elif p == 3:

                ydn1 = self.yd_func(i*h,self.y[i])
                ydn = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                ydn_2 = self.yd_func((i-3)*h,self.y[i-3])

                integral = h*((9.0*ydn1 + 19.0*ydn - 5.0*ydn_1 + ydn_2)/24.0)
        #if polinomial degree 4
        elif p == 4:

                ydn1 = self.yd_func(i*h,self.y[i])
                ydn  = self.yd_func((i-1)*h,self.y[i-1])
                ydn_1 = self.yd_func((i-2)*h,self.y[i-2])
                ydn_2 = self.yd_func((i-3)*h,self.y[i-3])
                ydn_3 = self.yd_func((i-4)*h,self.y[i-4])

                integral = ((251.0*ydn1)/720.0 + (646.0*ydn)/720.0 - (264.0*ydn_1)/720.0 + (106.0*ydn_2)/720.0 - (19.0*ydn_3)/720.0)*h


        y_tmp = self.y[i-1] + integral


        # Error very high! -> Decrease H
        if abs(y_tmp - self.y[i]) > self.episolon:
            print("Error very high! -> Decrease H")
            #self.h *= 0.9

            #self.h *= 0.8

        self.y[i] = y_tmp

        return self.y[i]

    # Prediction correction
    def __backward_diff__(self,p, h,i):


        y_tmp = 0


        if i < p:
            self.__runge_kutta__(h,i)
            return self.y[i]



        yi = symbols("yi")


        ydn1 = self.yd_expr.subs({self.y_symb:yi, self.x_symb:(i*h)})

        if p == 1:

            y_tmp = self.y[i-1] + h*ydn1

        elif p == 2:

            y_tmp = ((4.0/3.0)*self.y[i-1] - self.y[i-2]/3.0 + (2.0/3.0)*h*ydn1)

        elif p == 3:

            y_tmp = ((18.0/11.0)*self.y[i-1] - (9.0/11.0)*self.y[i-2] + (2.0/11.0)*self.y[i-3] + (6.0/11.0)*h*ydn1)

        elif p == 4:

            y_tmp = ((48.0/25.0)*self.y[i-1] - (36.0/25.0)*self.y[i-2] + (16.0/25.0)*self.y[i-3] - (3.0/25.0)*self.y[i-4] + (12.0/25.0)*h*ydn1)

        exp = y_tmp
        equ = Eq(exp,yi)

        self.y[i] = solve(equ,yi)[0]

        return self.y[i]




    def __select_method__(self,method="Euler"):
        
        method_func = self.__euler__

        if( method == "Euler" ):
            method_func = self.__euler__
        elif (method == "BackEuler"):
            method_func = self.__backward_euler__
        elif (method == "ImpEuler"):
            method_func = self.__improved_euler__
        elif (method == "RungeKutta"):
            method_func = self.__runge_kutta__
        elif (method == "Taylor"):

             # Partial derivative of yd w.r.t to x
            f_x = lambdify((self.x_symb,self.y_symb),diff(self.yd_expr,self.x_symb),"numpy")
            # Partial derivative of yd w.r.t to y
            f_y = lambdify((self.x_symb,self.y_symb),diff(self.yd_expr,self.y_symb),"numpy")

            method_func = partial(self.__3term_taylor_series__,f_x,f_y)

        elif (method == "Adams-Bashforth1"):

            method_func = partial(self.__adams_bashforth__,1)

        elif (method == "Adams-Bashforth2"):

            method_func = partial(self.__adams_bashforth__,2)
        elif (method == "Adams-Bashforth3"):

            method_func = partial(self.__adams_bashforth__,3)
        elif (method == "Adams-Bashforth4"):

            method_func = partial(self.__adams_bashforth__,4)

        elif (method == "Adams-Multon1"):

            method_func = partial(self.__adams_multon__,1)

        elif (method == "Adams-Multon2"):

            method_func = partial(self.__adams_multon__,2)
        elif (method == "Adams-Multon3"):

            method_func = partial(self.__adams_multon__,3)
        elif (method == "Adams-Multon4"):

            method_func = partial(self.__adams_multon__,4)

        elif (method == "Predictor-Corrector1"):

            method_func = partial(self.__predictor_corrector__,1)

        elif (method == "Predictor-Corrector2"):

            method_func = partial(self.__predictor_corrector__,2)
        elif (method == "Predictor-Corrector3"):

            method_func = partial(self.__predictor_corrector__,3)
        elif (method == "Predictor-Corrector4"):

            method_func = partial(self.__predictor_corrector__,4)

        elif (method == "BackDiff1"):

            method_func = partial(self.__backward_diff__,1)

        elif (method == "BackDiff2"):

            method_func = partial(self.__backward_diff__,2)
        elif (method == "BackDiff3"):

            method_func = partial(self.__backward_diff__,3)
        elif (method == "BackDiff4"):

            method_func = partial(self.__backward_diff__,4)




        return method_func

    def __solve__(self,x0,y0,h,n, method_func):

        self.h = h
        self.__initialize__(x0,y0,n)

        for i in range(1,n+1):

            self.x[i] = self.x[i-1] + self.h
            method_func(self.h,i)
            

            if self.phi_func != None:
                self.phi[i] = self.phi_func(self.h*i)
                self.error[i] = self.phi[i] - self.y[i]
                self.accerror[i] = abs(self.error[i]) + abs(self.error[i-1])

    def solve(self,x0,y0,h,n, method="Euler"):

        print( "------------ Solving [ yd =",self.yd_expr,"] for method:", method, "-------------")

        self.method = method
        method_func = self.__select_method__(method)


        self.__solve__(x0,y0,h,n,method_func)

    def plot(self,invert_yaxis = False ):
        plt.subplot(2, 1, 1)
        plt.title(self.method + " (h = " + str(self.h) +")")
        p1, = plt.plot(self.x, self.y, 'b', linewidth=1, label='y')
        p2, = plt.plot(self.x, self.phi, 'g', linewidth=1, label='phi(x)')
        plt.legend( [p1, p2], ['y', 'phi(x)'] )

        ax = plt.subplot(2, 1, 2)
        p3, = plt.plot(self.x,abs(self.error), 'r', linewidth=2, label='error')
        plt.legend( [p3], ['abs error'] )

        #verts = list(zip(self.x, abs(self.error)))

        #poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
        #ax.add_patch(poly)

        if invert_yaxis: 
            ax.invert_yaxis()
        plt.show()

    def plotExperiment(self, line_color = 'b', invert_yaxis = False ):
        #plt.clf()

        plt.subplot(2, 1, 1)
        plt.title(self.method)
        self.h = 0.1
        self.solve(0,1,self.h,int((1.0/self.h)),self.method)
        p1, = plt.plot(self.x, self.y, 'b', linewidth=1, label="y (h = " + str(0.1) +")")
        ax = plt.subplot(2, 1, 2)
        p51, = plt.plot(self.x,abs(self.error), 'b', linewidth=2, label='error')


        plt.subplot(2, 1, 1)
        self.h = 0.05
        self.solve(0,1,self.h,int((1.0/self.h)),self.method)
        p2, = plt.plot(self.x, self.y, 'c', linewidth=1, label="y (h = " + str(0.05) +")")
        ax = plt.subplot(2, 1, 2)
        p52, = plt.plot(self.x,abs(self.error), 'c', linewidth=2, label='error')


        plt.subplot(2, 1, 1)
        self.h = 0.025
        self.solve(0,1,self.h,int((1.0/self.h)),self.method)
        p3, = plt.plot(self.x, self.y, 'm', linewidth=1, label="y (h = " + str(0.025) +")")
        ax = plt.subplot(2, 1, 2)
        p53, = plt.plot(self.x,abs(self.error), 'm', linewidth=2, label='error')
        plt.legend( [p51, p52, p53], ['abs error (h = 0.1)','abs error (h = 0.05)', 'abs error (h = 0.025)'], loc = 0 )


        plt.subplot(2, 1, 1)
        p4, = plt.plot(self.x, self.phi, 'g', linewidth=2, label='phi(x)')
        plt.legend( [p1, p2,p3, p4], ["y (h = 0.1)","y (h = 0.05)","y (h = 0.025)", 'phi(x)'], loc=0 )


        #verts = list(zip(self.x, abs(self.error)))

        #poly = Polygon(verts, facecolor='0.9', edgecolor='0.5')
        #ax.add_patch(poly)

        if invert_yaxis: 
            ax.invert_yaxis()
        #plt.show()
        plt.savefig("plots/"+self.method+ ".png")



def readInput(filename):

    yd_str = None
    yd_expression = None
    phi_str = None
    phi_expression = None
    y0 = 0
    n = 10
    h = 0.1
    episolon = 0.1
    method = "Euler"

    try:

        input_file = open(filename, 'r')
        
        lines = input_file.readlines()

        valid_lines = [line.replace(' ','').replace('\n','').replace('\r','') for line in lines if ((not line.startswith("#")) and len(line.replace(" ",'')) > 1)]

   
        # y' = g(x) - p(x)*y
        yd_str = valid_lines[0].split("=")[1]

        yd_expression = yd_str

        # phi(x): analytical solution for error comparisson
        phi_str = valid_lines[1].split("=")[1]
        phi_expression = phi_str if len(phi_str) > 1 else None



        y0_str = valid_lines[2]
        h_str = valid_lines[3]
        n_str = valid_lines[4]
        episolon_str = valid_lines[5]

        y0 = float(y0_str.split("=")[1])
        h = float(h_str.split("=")[1])
        n = float(n_str.split("=")[1])

        episolon = float(episolon_str.split("=")[1])

        method = valid_lines[6]
    except:
        print( "WRONG/NON-EXISTENT FILE NAME OR MALFORMED INPUT FILE")
        sys.exit(1)

    return (yd_expression,phi_expression,y0,n,h,episolon,method)


def main(argv=None):

    yd_expression = None
    phi_expression = None
    y0 = 0
    n = 0
    h = 0.1
    method = "Euler"

    yd_expression,phi_expression,y0,n,h,episolon,method = readInput(argv[1] if len(argv)>1 else "inputFile.txt")
    
    ds = DSolver(yd_expression, phi_expression,episolon)

    ds.solve(0,y0,h,int((1.0/h)+1),method)

    print( "Y: ", ds.y)
    print( "Phi: ", ds.phi)
    print( "Error: ", ds.error)
    print( "Acumulated Error: ", sum(abs(ds.error)))

    ds.plot()


def mainExperiment(argv=None):

    starting_methods = ["Euler", "BackEuler", "ImpEuler", "RungeKutta", "Taylor"]
    multistep_methods = ["Adams-Bashforth", "Adams-Multon","Predictor-Corrector","BackDiff"]
    colors = ['b','m','c']
    hs = [0.1, 0.05, 0.025]
    yd_expression = None
    phi_expression = None
    y0 = 0
    n = 0
    h = 0.1
    method = "Euler"

    yd_expression,phi_expression,y0,n,h,episolon,method = readInput(argv[1] if len(argv)>1 else "inputFile.txt")
    
    ds = DSolver(yd_expression, phi_expression,episolon)

    out_file = open("experiment_report.txt",'w')

    
    for method in starting_methods:
        plt.clf()
        out_file.write("----------------------------------------------------------------\n")
        out_file.write("METHOD: " + method + "\n")


        x = y = phi = error = accerror = []
        for j in range(0,len(hs)):
            h = hs[j]
            ds.solve(0,y0,h,int((1.0/h)),method)
            out_file.write("H: " + str(h) + " N: " + str(int((1.0/h))) + "\n")
            out_file.write("X: " + str(ds.x) + "\n")
            out_file.write("Y: " + str(ds.y) + "\n")
            out_file.write("Phi: " + str(ds.phi)+ "\n")
            out_file.write("Error: " + str(ds.error)+ "\n")
            out_file.write("Acumulated Error: " + str(sum(abs(ds.error))) + "\n")

            


        ds.plotExperiment(colors[j])

    for k in range(0,len(multistep_methods)):
        for i in range(1,5):
            plt.clf()
            method = multistep_methods[k] + str(i)
            out_file.write("METHOD: " + method + "\n")
            for j in range(0,len(hs)):
                h = hs[j]
                ds.solve(0,y0,h,int((1.0/h)),method)
                out_file.write("H: " + str(h) + " N: " + str(int((1.0/h))) + "\n")
                out_file.write("X: " + str(ds.x) + "\n")
                out_file.write("Y: " + str(ds.y) + "\n")
                out_file.write("Phi: " + str(ds.phi)+ "\n")
                out_file.write("Error: " + str(ds.error)+ "\n")
                out_file.write("Acumulated Error: " + str(sum(abs(ds.error)))+ "\n")


            ds.plotExperiment(colors[j])
            #out_file.write('\n')


        #out_file.write('\n')
    #out_file.close()




    #print "Y: ", ds.y
    #print "Phi: ", ds.phi
    #print "Error: ", ds.error
    #print "Acumulated Error: ", sum(abs(ds.error))

    #ds.plot()

    out_file.close()

if __name__ == "__main__":
    sys.exit(main(sys.argv))