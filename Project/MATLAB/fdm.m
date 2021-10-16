%%% Solving of Partial Differential Equation using
%%% Finite Difference Method

clc
clear
close all
n0 = 1;
a = 1;
tmin = 0 ;
tmax = 10;
nt = 100;
dt = 1/nt;
t = tmin:dt:tmax;
n = n0*exp(-a*t);
plot(t,n,'Linewidth',2)

grid on
hold on

Nt = floor(tmax-tmin)/dt;
N = zeros(1,Nt);
T = zeros(1,Nt);

N(1) = 1;
for i =1:Nt 
    T(i) = i*dt;
    N(i+1) = N(i)*(1-a*dt); 
end
TT = [0,T]
plot(TT,N,'--','Linewidth',2)

%% Error
figure()
E = n-N;
plot(E,'Linewidth',2)

