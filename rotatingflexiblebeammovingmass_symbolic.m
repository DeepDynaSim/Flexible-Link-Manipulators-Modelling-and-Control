% Corrected and fully optimized MATLAB code including the definition of Kf

% Define variables - Optimization does not change variable definitions
syms x L rho E A Iz th thdot thddot q1 q1dot q1ddot q2 q2dot q2ddot m u udot uddot g real;

% Precompute expressions to reduce redundant computations
x_over_L = x/L;
x_over_L2 = x_over_L^2;
x_over_L3 = x_over_L^3;
u_over_L2 = (u/L)^2;
u_over_L3 = (u/L)^3;
u_over_L4 = (u/L)^4;
u_over_L5 = (u/L)^5;
u_over_L6 = (u/L)^6;

% Optimized shape functions for beam element
N1 = 3*x_over_L2 - 2*x_over_L3;
N2 = x_over_L3 - x_over_L2;
N = [N1, N2]; % Form shape function matrix efficiently

% Define stiffness matrix Kf for the beam
Kf = E*Iz/L^3 * [12, -6*L; -6*L, 4*L^2];

% Compute elastic mass matrix and mass coupling term once, avoiding redundancy
Mf = rho*A*int(N'*N, x, 0, L);
Mrf = rho*A*int(N*x, x, 0, L);

% Vector of elastic coordinates and velocities
q = [q1; q2];
qdot = [q1dot; q2dot];

% Precompute Nu for particle dynamics to avoid redundant subs call
Nu = subs(N, x, u);

% Kinetic Energy (KE) computations optimized by precomputing reused expressions
T_beam_part = thdot^2*q'*Mf*q + qdot'*Mf*qdot + 2*thdot*Mrf*qdot;
T_particle_part = m*(udot^2 - 2*udot*thdot*Nu*q + thdot^2*q'*Nu'*Nu*q + qdot'*Nu'*Nu*qdot + 2*thdot*u*Nu*qdot + thdot^2*u^2);
T = 1/2*(T_beam_part + T_particle_part);

% Potential Energy (PE) computations optimized
Md = rho*A*int(N, x, 0, L);
V_beam = 1/2*q'*Kf*q + rho*A*g*L^2*sin(th)/2 + g*cos(th)*Md*q;
V_particle = m*g*(u*sin(th) + cos(th)*Nu*q);
V = V_beam + V_particle;

% Lagrangian of the system
Lag = T - V;

% Equations of Motion (EoM) optimization by grouping similar derivative operations
S = [diff(Lag, q1dot), diff(Lag, q2dot), diff(Lag, udot), diff(Lag, thdot)];
E = [diff(S, q1) * q1dot + diff(S, q2) * q2dot + diff(S, q1dot) * q1ddot + diff(S, q2dot) * q2ddot + diff(S, u) * udot + diff(S, udot) * uddot + diff(S, th) * thdot + diff(S, thdot) * thddot];
D = [diff(Lag, q1), diff(Lag, q2), diff(Lag, u), diff(Lag, th)];

% Form LHS of EoM for q1, q2, u, and th optimized
LHS = simplify(E - D);

% Displaying the optimized equations of motion
disp('LHS of Equation of Motion for q1:');
disp(LHS(1));
disp('LHS of Equation of Motion for q2:');
disp(LHS(2));
disp('LHS of Equation of Motion for u:');
disp(LHS(3));
disp('LHS of Equation of Motion for th:');
disp(LHS(4));