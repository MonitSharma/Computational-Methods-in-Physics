%%% Analytical Solutiopn of 2D Laplace Equation


clc
clear
close all

a = 3*pi;
b = pi;

nx = 100;
ny = 50;
dx = a/nx;
dy = b/ny;


x = 0:dx:a;
y = 0:dy:b;

[X,Y] = meshgrid(x,y)
u = sin(X)./sin(a)*sinh(Y)./sinh(b);
maxu = max(u,[],'all')
nu = u/maxu;
surf(u)
colormap(jet(255))

shading interp
%figure()
%contour(nu)