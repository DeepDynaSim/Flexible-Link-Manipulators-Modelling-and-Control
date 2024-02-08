% Ensure you have the Symbolic Math Toolbox installed and available

syms q1(t) q2(t) % Defining symbolic functions of time for the motor angle and the second mode

syms I_h I_b rho omega_2 gamma_2 tau % Defining symbolic constants

% Define the derivatives of q1 and q2 for use in the Lagrangian
dq1 = diff(q1, t);
dq2 = diff(q2, t);

% Define the Lagrangian L for the second mode only
L = 1/2 * (I_h + I_b) * dq1^2 + rho/2 * dq2^2 + rho * dq1 * dq2 * gamma_2 - rho/2 * omega_2^2 * q2^2;

% Compute the Euler-Lagrange equations for q1 and q2

% For q1
dL_dq1 = diff(L, q1); % Correctly use diff here since functionalDerivative is not needed
dL_ddq1 = diff(L, diff(q1, t)); % Correct approach for derivatives
EL_q1 = diff(dL_ddq1, t) - dL_dq1 == tau;

% For q2
dL_dq2 = diff(L, q2); % Similar correction as for q1
dL_ddq2 = diff(L, diff(q2, t)); % Correct approach for derivatives
EL_q2 = diff(dL_ddq2, t) - dL_dq2 == tau; 

% Display the Euler-Lagrange equations for q1 and q2
disp('Euler-Lagrange Equation for q1:');
disp(EL_q1);
disp('Euler-Lagrange Equation for q2:');
disp(EL_q2);