% Define symbolic variables and equation

syms phi(x) omega 

%EI rho A L P

% Given properties

a = 6.35e-3; % meters
A = a^2; % Cross-sectional area
I = a^4 / 12; % Moment of inertia
E = 69e9; % Pa
rho = 2700; % kg/m^3
L = 1; % meters
P = -10; % Newtons
EI=E*I;

% Beam equation with omega squared term on the right-hand side
D_eqn = diff(phi, x, 4) == omega^2 * phi*(rho*A/EI);

% Attempt to define boundary conditions symbolically
BCs = [phi(0) == 0, ...
       subs(diff(phi, x), x, 0) == 0, ...
       subs(diff(phi, x, 2), x, L) == P*L/EI, ...
       subs(diff(phi, x, 3), x, L) == 0];

% Solve the differential equation
phiSol(x) = dsolve(D_eqn, BCs);

% Display the solution
disp(phiSol);

% First mode approximated root for clamped-free beam
k1 = 1.87510407;
% Calculate omega_1 using the given relationship
omega_1 = k1^2 * sqrt(E*I/(rho*A));
phiSol1=subs(phiSol,omega,omega_1);
syms x_sym
phiSol1=subs(phiSol1,x,x_sym);

phiSol1=subs(phiSol1,x_sym,x);
% Convert symbolic solution to function handle for plotting
phiSol1_func = matlabFunction(phiSol1);

% Generate x values from 0 to L for plotting
x_vals = linspace(0, L, 100);

% Evaluate phiSol1 over the range of x values
phiSol1_vals = phiSol1_func(x_vals);

% Plot the solution
figure(1);
plot(x_vals, phiSol1_vals, 'r-', 'LineWidth', 2);
xlabel('x (meters)');
ylabel('\phi(x)');
title('Beam Mode Shape for the First Approximated Root');
grid on;