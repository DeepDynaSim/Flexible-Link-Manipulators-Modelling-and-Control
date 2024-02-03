% MATLAB Script for Analyzing Mode Shapes of a Clamped-Free Beam with Modified BC

% Define symbolic variables and parameters
syms w(x) C1 C2 C3 C4 EI L M0
assume(EI > 0);
assume(L > 0);
assume(M0, 'real'); % Assuming M0 represents an applied moment at the free end due to step loading

% Define the differential equation according to Euler-Bernoulli beam theory
diffEq = EI*diff(w, x, 4) == 0;

% Solve the differential equation for w(x)
wSol(x) = dsolve(diffEq);

% Apply clamped-free boundary conditions with an assumed moment M0 at the free end:
% w(0) = 0, w'(0) = 0 (clamped conditions)
% EI*w''(L) = M0, EI*w'''(L) = 0 (free end with applied moment M0 and zero shear)
BCs = [subs(wSol, x, 0) == 0, ...
       subs(diff(wSol, x), x, 0) == 0, ...
       subs(EI*diff(wSol, x, 2), x, L) == M0, ...
       subs(EI*diff(wSol, x, 3), x, L) == 0];

% Solve for the constants C1, C2, C3, C4
constants = solve(BCs, [C1, C2, C3, C4]);

% Substitute constants back into the general solution
modeShape = subs(wSol, constants);

% Display the mode shape equation
disp('Mode Shape Equation with Modified Boundary Conditions:');
disp(modeShape);

% Optional: Plot the first mode shape for a specific applied moment M0
fplot(subs(modeShape, {EI, L, M0}, {69e9 * (6.35e-3)^4 / 12, 1, 10}), [0, 1]); % Example values EI, L, M0 for plotting
title('Mode Shape of Clamped-Free Beam with Applied Moment at Free End');
xlabel('x (m)');
ylabel('Deflection w(x)');
grid on;