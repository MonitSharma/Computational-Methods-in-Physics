%%% Numerical Solution of 1D Wave equation
%%% using Finite difference technique

clc  
clear
close all





nx = 100;
nt = 300;
a = 0;
b = 1;
t0 = 0;
tf = 2;
T0 = 40;
rho0 = 0.01;
alpha = 0.5;

dx = (b-a)/(nx-1);
dt = (tf - t0)/(nt - 1);

x = a:dx:b;
t = t0:dt:tf;

s = dt^2/ dx^2

%% Analytical SOlution


UA = zeros(nx,nt);
for i=1:nx
    for j = 1:nt
        UA(i,j) = sin(pi*x(i))* (cos(pi*t(j)) + sin(pi*t(j))/pi);
    end
end


figure()
countourf(UA,200,'linecolor', 'non')
title('Analytical Solution')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar




%% Numerical solution

UN = zeros(nx,nt);


UN(:,1) = sin(pi.*x);
UN(:,2) = sin(pi.*x)*(1 + dt);


for j = 2:nt-1
    for i = 2:nx-1
        UN (i, j+1) = 2*UN(i,j) + UN(i,j-1) + (dt^2/dx^2) * (rho0/T0)(UN(i+1,j) - 2*y(i,j) + UN(i-1,j)+ alpha*T0/rho0(UN(i+1,j)- UN(i,j)dx));
    end
end


figure()
countourf(UN,200,'linecolor', 'non')
title('Numerical Solution')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar

%% Error

E = abs(UA-UN);

figure()
countourf(E,200,'linecolor', 'non')
title('Error')
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar
