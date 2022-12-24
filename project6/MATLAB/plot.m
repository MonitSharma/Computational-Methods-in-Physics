clc
clear
close all

nx = 100;
nt = 500;
a = 0;
b = 1;
t0 = 0;
tf = 50;
alpha = 0.5;
rho = 0.01;
ten = 40.0;
c = rho/ten;
l= 1;



dx = (b-a)/(nx-1);

dt = (tf-t0)/(nt-1);

x = a:dx:b;
t = t0:dt:tf;

s = dt^2/dx^2;

%% Analytical Solution
UA = zeros(nx,nt);
for i =1:nx
    for j = 1:nt
        UA(i,j) = sin(i*pi*x(i)/l)*sin(pi*x(i)/l)*cos(c*pi*t(j)/l);
    end
end
figure()
contourf(UA,200,'linecolor','non');
title('Analytical Solution');
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar



%% Numerical Solution

UN = zeros(nx,nt);

% Boundary Conditions
for j = 1:nt-1
    for i = 1:nx-1
        if x(i) <= 0.5
            UN(i,1) = 1.25*x(i)*l;
        else 
            UN(i,1) = 5*(1-x(i)/l);
        end
    end
end

%UN(:,1) = 1.25*x/l;
%UN(:,1) = 5*(1-x/l);
%UN(:,2) = sin(pi*x)*(1+dt);
UN(1,:) = 0;
UN(l,:) = 0;


for j = 2:nt-1
    for i = 2:nx-1
        UN(i,j+1) = c*s*(UN(i+1,j) + UN(i-1,j)-2*UN(i,j)) + alpha*c*(UN(i+1,j) - UN(i,j)) - UN(i,j-1) + 2*UN(i,j);
       
        %UN(i,j+1) = s*(UN(i-1,j) -2*UN(i,j) + UN(i+1,j))+ 2*UN(i,j) - UN(i,j-1);
    end
end




figure()
contourf(UN,200,'linecolor','non');
title('Numerical Solution');
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar


%% Error

E = abs(UA - UN);
figure()
contourf(E,200,'linecolor','non');
title('Error');
xlabel('t')
ylabel('x')
colormap(jet(256))
colorbar





