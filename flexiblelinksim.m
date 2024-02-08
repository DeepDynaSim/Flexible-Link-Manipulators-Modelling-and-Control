function flexiblelinksim
    % Constants
    I_h = 0.3; % Moment of inertia
    rho = 2700; % Density (or some specific parameter related to the system)
    omega = 32.5819; % Natural frequency of q2
    gamma = 0.1521*10^-10; % Damping factor evaluated from some specific conditions

    % Initial conditions for q1, q2, and their derivatives
    q0 = [0; 0; 0; 0]; % [q1; q1dot; q2; q2dot]

    % Time span
    tspan = [0, 1]; % 1 seconds simulation

    % Options for ODE solver to control the step size for more stability
    options = odeset('MaxStep', 1e-3);

    % ODE solver
    [t, x] = ode45(@(t,x) stateSpace(t, x, I_h, rho, omega, gamma), tspan, q0, options);

    % Plot results for positions and velocities
    figure;
    subplot(2,2,1);
    plot(t, x(:,1));
    title('q1 vs. Time');
    xlabel('Time (s)');
    ylabel('Position (rad)');

    subplot(2,2,2);
    plot(t, x(:,2));
    title('q1dot vs. Time');
    xlabel('Time (s)');
    ylabel('Velocity (rad/s)');

    subplot(2,2,3);
    plot(t, x(:,3));
    title('q2 vs. Time');
    xlabel('Time (s)');
    ylabel('Position (rad)');

    subplot(2,2,4);
    plot(t, x(:,4));
    title('q2dot vs. Time');
    xlabel('Time (s)');
    ylabel('Velocity (rad/s)');
end

function dxdt = stateSpace(t, x, I_h, rho, omega, gamma)
    q1 = x(1);
    q1dot = x(2);
    q2 = x(3);
    q2dot = x(4);

    % Smooth step response for tau using a sigmoid function for smoother transition
    tau = 0.25 / (1 + exp(-10 * (t - 0.5))); % Smoothly transitioning tau around t=1 second

    % Equations for q1 and q2
    q1ddot = (tau + rho*omega^2*gamma*q2) / I_h;
    q2ddot = (-tau*gamma/I_h - q2^2*omega^2*(1 + (rho*gamma^2)/I_h));

    dxdt = [q1dot; q1ddot; q2dot; q2ddot];
end