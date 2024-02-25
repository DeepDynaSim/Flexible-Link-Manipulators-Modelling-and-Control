% Comprehensive MATLAB Script for Acrobot System Dynamics and Zero-Dynamics Analysis

%% Define the symbolic variables 
syms theta1(t) theta2(t) tau(t) real 
syms m1 m2 l1 l2 g
%% Define the derivatives
dtheta1=diff(theta1,t);
dtheta2=diff(theta2,t);
%% Kinetic and Potential Energy
T = 0.5 * m1 * (l1^2) * dtheta1^2 + ...
    0.5 * m2 * ((l1 * dtheta1)^2 + (l2 * dtheta2)^2 + 2 * l1 * l2 * dtheta1 * dtheta2 * cos(theta2));
V = -m1 * g * l1 * cos(theta1) - m2 * g * (l1 * cos(theta1) + l2 * cos(theta1 + theta2));

%% Lagrangian
L = T - V;

%% Equations of Motion (Euler-Lagrange Equations)
%eq1 = diff2(diff2(L, dtheta1), t) - diff2(L, theta1) == 0;
%eq2 = diff2(diff2(L, dtheta2), t) - diff2(L, theta2) == tau;

%% Define system parameters

m1 = 8; % Mass of the first link in kg
m2 = 8; % Mass of the second link in kg
l1 = 0.5; % Length of the first link in meters
l2 = 1; % Length of the second link in meters
g = 10; % Acceleration due to gravity in m/s^2

syms ddth1 ddth2 dth1 dth2 theta1 theta2 tau

eq1=(m2*(2*l1^2*ddth1 - 2*l1*l2*sin(theta2)*dth2^2 + 2*l1*l2*cos(theta2)*ddth2))/2 + g*m2*(l1*sin(theta1) + l2*sin(theta1 + theta2)) + l1^2*m1*ddth1 + g*l1*m1*sin(theta1) == 0;

eq2=(m2*(2*l2^2*ddth2 + 2*l1*l2*cos(theta2)*ddth1 - 2*l1*l2*sin(theta2)*dth1*dth2))/2 + g*l2*m2*sin(theta1 + theta2) == tau;
 

%% Solve for ddtheta1 and ddtheta2
[ddtheta1_sol, ddtheta2_sol] = solve([eq1, eq2], [ddth1, ddth2]);

%% Define the ODE function for full system dynamics

syms t
syms y [4 1] real; % Define y as a column vector of 4 symbolic variables

% Extract components from y for clarity in defining the ODE system
theta1 = y(1);
theta2 = y(2);
dth1 = y(3);
dth2 = y(4);

% Define the step response for tau
% tau = 1 Nm for t > 1 second, and tau = 0 for t <= 1 second
tau = heaviside(t - 1); % heaviside function steps from 0 to 1 at t = 1

% Define ddth1, ddth2 using the extracted components and the step response for tau
ddth1 = (- 8*sin(theta2)*dth2^2 + 4*dth1*cos(theta2)*sin(theta2)*dth2 + 160*sin(theta1 + theta2) + 160*sin(theta1) + tau*cos(theta2) - 80*sin(theta1 + theta2)*cos(theta2))/(4*(cos(theta2)^2 - 2));
ddth2 = -(- 4*cos(theta2)*sin(theta2)*dth2^2 + 4*dth1*sin(theta2)*dth2 + tau - 80*sin(theta1 + theta2) + 80*cos(theta2)*sin(theta1) + 80*sin(theta1 + theta2)*cos(theta2))/(4*(cos(theta2)^2 - 2));

% Convert the system of equations to a MATLAB function handle
% Note: Since tau is now a function of t, we don't pass tau as a separate variable
odeFun = matlabFunction([dth1; dth2; ddth1; ddth2], 'Vars', {t, y});

%% Initial Conditions and Time Span
initial_conditions = [pi/4; 0; 0; 0]; % [theta1; theta2; dtheta1; dtheta2]
t_span = [0, 10]; % Time span for the simulation

%% Simulate Full System Dynamics
% Since tau is defined inside the ODE, we do not need to pass it separately
[t, Y] = ode45(@(t,y) odeFun(t,y), t_span, initial_conditions);

%% Plot Full System Dynamics
figure;
subplot(2,1,1);
plot(t, Y(:,1), 'DisplayName', 'Theta1');
hold on;
plot(t, Y(:,2), 'DisplayName', 'Theta2');
title('Full System Dynamics of Theta1 and Theta2 over Time');
xlabel('Time (s)');
ylabel('Angle (rad)');
legend;
grid on;

%% Linearization

% Define Symbolic Variables
syms theta1 theta2 dth1 dth2 tau real

% System Equations
ddth1 = (- 8*sin(theta2)*dth2^2 + 4*dth1*cos(theta2)*sin(theta2)*dth2 + 160*sin(theta1 + theta2) + 160*sin(theta1) + tau*cos(theta2) - 80*sin(theta1 + theta2)*cos(theta2))/(4*(cos(theta2)^2 - 2));
ddth2 = -(- 4*cos(theta2)*sin(theta2)*dth2^2 + 4*dth1*sin(theta2)*dth2 + tau - 80*sin(theta1 + theta2) + 80*cos(theta2)*sin(theta1) + 80*sin(theta1 + theta2)*cos(theta2))/(4*(cos(theta2)^2 - 2));

% System Dynamics Vector
f = [dth1; dth2; ddth1; ddth2];

% State and Input Vectors
y = [theta1; theta2; dth1; dth2];
u = tau;

% Calculate Jacobian Matrices A and B
A = jacobian(f, y);
B = jacobian(f, u);

% Evaluate Jacobian Matrices at the Equilibrium Point
A_lin = subs(A, {theta1, theta2, dth1, dth2, tau}, {0, 0, 0, 0, 0});
B_lin = subs(B, {theta1, theta2, dth1, dth2, tau}, {0, 0, 0, 0, 0});

% Assuming the linearized matrices A_lin and B_lin are obtained from the previous steps

% Use the actual linearized system matrix A_lin and input matrix B_lin
A = double(A_lin); % Convert from symbolic to numeric matrix for ss function
B = double(B_lin);

% Outputs are theta1 and theta2
C = [1 0 0 0; 0 1 0 0]; % Output matrix
D = [0; 0]; % Direct transmission matrix

% Create a state-space model with the linearized system
sys = ss(A, B, C, D);

% Time vector for simulation
t = linspace(0, 10, 1000); % More points for smoother step transition

% Define tau(t) as a step function with step at 1 second and amplitude of 1 Nm
tau = ones(size(t)); % Initialize tau with amplitude of 1 Nm
tau(t < 1) = 0; % Set tau to 0 for t < 1 second

% Initial condition for the state vector
initialCondition = zeros(size(A, 1), 1); % Assuming starting from equilibrium

% Simulate the response of the linearized system to the step input in tau
[y, t, x] = lsim(sys, tau, t, initialCondition);

%% Plot the Response of the Linearized System to Step Input in Tau
subplot(2,1,2);
plot(t, x);
title('Response of Linearized System to Step Input in Tau');
xlabel('Time (s)');
ylabel('State Response');
legend('Theta1', 'Theta2', 'dTheta1', 'dTheta2');
grid on;