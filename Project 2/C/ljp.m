clc
clear
close all

Eb = 10.4;     %Binding energy in meV
r_o = 3.84;     % equilibrium distance in A
syms A B;       % A & B are LJ constants
[A,B] = solve((A^2)/(4*B) == Eb, (2*B/A)^(1/6)==r_o,A,B);  % to solve A & B

A = double(A);
B = double(B);


r = 0:0.1:8;     % interatomic distance
phi_r = -A./r.^6+B./r.^12  ; % LJ Potential Function

phi_o = A^2/(4*B) ;     % potential at equilibrium distance


fr = 6*A./r.^7-12*B./r.^13;    % force attraction

 plot(r/r_o, phi_r/phi_o,'LineWidth',1)
 xlabel('r/r_o')
 ylabel('\phi(r)/\phi_o')
 ylim([-1.5 2])


 hold on
 plot(r/r_o, fr/phi_o,'--','LineWidth',1)
 legend('L-J potential', 'force')
 ylim([-1.5 2])


 grid on
