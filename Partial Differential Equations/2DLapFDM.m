%%% Numerical Solution of 2D Laplace Equation using
%%% Finite Difference Method for Iterative Techique

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
u = zeros(ny+1, nx+1);
for i = 1:ny+1
    for j = 1:nx+1
        u (i,j) = sin(x(j))./sin(a).*sinh(y(i))./sinh(b);
    end
end

maxu = max(u,[], 'all');
nu = u/maxu;

contourf(nu,200,'linecolor','none');
title('Analytical')
xlabel('x')
ylabel('y')
colormap(jet(255));

colorbar
caxis([-1,1])



%%% Numerical  functions
Co = 1/(2*(dx^2+dy^2));
% initial guess
U = zeros(ny+1, nx+1);

% Boundary Condition

U (1,:) = 0;
U (ny+1,:) = sin(x)/sin(a);
U (:,1) = 0;

U (:, nx+1) = sinh(y)/sinh(b);
iter_number = 10;
for k =1:iter_number

for i = 2:ny
    for j = 2:nx
        U(i,j) = Co*(dx^2(U(i+1,j)+ U(i-1,j))+dy^2*(U(i,j+1)- U(i,j-1)));
    end
end
end

maxU = max(U,[], 'all');
nU = U/maxU;
figure()
contourf(nU,200,'linecolor','none');
title('Numerical')
xlabel('x')
ylabel('y')
colormap(jet(255));

colorbar
caxis([-1,1])

%%% Error    calculation
E = nu-nU;
figure()
contourf(E,200,'linecolor','none');
title(Error')
xlabel('x')
ylabel('y')
colormap(jet(255));