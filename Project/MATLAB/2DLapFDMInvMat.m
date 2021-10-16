%%%% Numerical Solution of 2D Laplace equation using
%%% FDM and Inverse Matrix Technique



clc
clear 
close all

nx = 100;
ny = 50;
k = nx*ny;
a = 3*pi;
b = pi;
dx = a/(nx+1);
dy = b/(ny+1);
x = 0:dx:a;
y = 0:dy:b;


alpha = (dx/dy)^2;
B = 2*(1+alpha);

array1 = ones(1,k);
array2 = ones(1,k-1);
array3 = ones(1,k-nx);


Adiag = diag(array1)*-B;
Aupdiag1 = diag (array2,1);
Aupdiag2 = diag (array3,nx)* alpha;


if nx == ny

for i = 1:nx-1
    Aupdiag1(i*nx, i*nx+1)= 0;
end
end


if nx < ny

    for i = 1:ny-1
        Aupdiag1(i*nx, i*nx+1)= 0;
    end
end


if nx > ny

    for i = 1:nx-(nx-ny+1)
        Aupdiag1(i*nx, i*nx+1)= 0;
    end
end



Aup = Aupdiag1 + Aupdiag2;

Adown = Aup.';

A1 = Aup + Adown;

A = Adiag + A1;

%%

U = zeros(ny+2, nx+2);

%%% BC

U(1,:) = 0;
U(ny+2,:) = sin(x)/sin(a);
U(:,1) = 0;
U(:, nx+2)= sinh(y)/sinh(b);
U;

F = zeros(k,1);
c = 1;
for i = 1:ny
    for j = 1:nx
        if j ==1
            F(c) = F(c) - U(i+1,j);
        end
        if i ==1
            F(c) = F(c) - alpha*U(i,j+1);
        end
        if j ==nx
            F(c) = F(c) - U(i+1,j+2);
        end

        if i ==ny
            F(c) = F(c) - alpha*U(i+2,j+1);

        end
        c = c+1;
    end
end



u = linsolve(A,F);

ur = zeros(ny,nx);
p = 1;
for i = 1:ny
    for j = 1:nx
        ur (i,j) = u(p);
    p = p+1;
    end
end




for i = 1:ny
    for j = 1:nx
        U (i+1,j+1) = ur(i,j);
    
    end
end


UN = U;
maxUN = max(UN,[],'all');
UN = UN/maxUN;
figure()

contourf(UN,200,'linecolor','non')
title('Numerical')
xlabel('x')
ylabel('y')
colormap(jet(256))
caxis([-1,1])
savefig('Numerical.fig')

%%% analytical Solution
UA = zeros(ny+2,nx+2);
for i = 1:ny+2
    for j = 1:nx+2
        UA(i,j) = sin(x(j))./sin(a)*sinh(y(i))./sinh(b);
    end
end

maxUA = max(UA,[],'all');
UA = UA/maxUA;

figure()
contourf(UA,200,'linecolor','non')
title('Analytical')
xlabel('x')
ylabel('y')
colormap(jet(256))
caxis([-1,1])
savefig('Analytical.fig')


%%% Error
 E = UA - UN;
 figure()
 contourf(E,200,'linecolor','non')
 title('Error ')
 xlabel('x')
 ylabel('y')

 colormap(jet(256))
 colorbar
 caxis([-1,1])
 savefig('Error.fig')

