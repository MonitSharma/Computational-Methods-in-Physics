#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define min(a, b) ((a) < (b) ? (a) : (b)) 

double f (double t, double y);

int main ()
{
	double x0,y0,k1,k2,k3,k4,k5,k6,w,w1,w2,y,x,xn,y1,y2,h,R,epsilon,delta, g;
	int i=0;
	
	//boundary conditions 
	epsilon = 0.00001;
	x0=1;
	y0=2; 
	xn=3; 
	h=0.2; 
	x=x0;
	y=y0;
	
	
	printf("Step %d: x=%f\t y=%f\n", i, x, y);
	while (x<=xn)
	{
		h = min(h, 3-x);
			
		k1 = h*f(x,y);
		k2 = h*f(x+h/4, y+k1/4);	
		k3 = h*f(x+3*h/8,  y+3*k1/32+9*k2/32);		
		k4 = h*f(x+12*h/13, y+1932*k1/2197-7200*k2/2197+7296*k3/2197);
		k5 = h*f(x+h, y+439*k1/216-8*k2+3680*k3/513-845*k4/4104);
		k6 = h*f(x+h/2, y-8*k1/27+2*k2-3544*k3/2565+1859*k4/4104-11*k5/40);
		
		y1 = y + (25*k1/216+1408*k3/2565+2197*k4/4104-k5/5);
		y2 = y + (16*k1/135+6656*k3/12825+28561*k4/56430-9*k5/50+2*k6/55);
		
		R = ((fabs(y2-y1)))/h;
		delta = 0.84*(pow((epsilon/R), 0.25));
		
		if (R<=epsilon)
		{
			x = x+h;
			y = y1;
			i = i+1;
			printf("Step %d: x=%f\t y=%f\n", i, x, y);
			h = delta*h;
		}
		else
		{
			h = delta*h;
            x = x+h;	
		}
	}
}

double f(double x, double y)
{
	double v;
	v = pow(x,-2)*(x*y - pow(y,2));
	return v;
}