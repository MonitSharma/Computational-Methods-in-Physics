%%% Numerical Solution of 1D Heat Equation using Finite 
%%% Differemce method


clc
clear
close all

nx = 10;
nt = 300;
a = 0;
b = 1;
t0 = 0;
tf = 0.2;
dx = (b-a)/(nx-1);
dt = (tf-t0)/(nt-1);
x = a:dx:b;
t = t0:dt:tf;

s = dt/dx^2;

%% Analytical Solution
UA = zeros(nx,nt);

for j = 1:nt
    for i = 1:nx
        UA(i,j) = sin(pi*x(i))*exp(-pi^2*t(j));
    end
end

figure()

contourf(UA,200,'linecolor','non')
xlabel('x')
ylabel('t')
title('Analytical Solution')
colormap(jet(256))

colorbar
caxis([0,1])



%%% Numerical Solution

% initial cond


UN = zeros(nx,nt);
UN(:,1) = sin(pi*x);
for j = 1:nt-1
    for i = 2:nx-1
        UN(i,j+1) = s*UN(i-1,j) + (1-2*s)*UN(i,j) + s*UN(i+1,j);
    end
end

figure()

contourf(UN,200,'linecolor','non')
xlabel('x')
ylabel('t')
title('Numerical Solution')
colormap(jet(256))

colorbar
caxis([0,1])


% Error
E = abs(UA-UN)

figure()

contourf(E,200,'linecolor','non')
xlabel('x')
ylabel('t')
title('Error')
colormap(jet(256))

colorbar
caxis([0,1])
