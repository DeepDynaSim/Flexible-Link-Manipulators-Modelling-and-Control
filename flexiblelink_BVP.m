function flexiblelink_BVP
    % Define the beam properties
    L = 1; % Length of the beam in meters
    EI = 69e9 * (6.35e-3)^4 / 12; % Young's modulus times moment of inertia
    M0 = 10; % Applied moment at the free end (N*m), matching the analytical example

    % Define the differential equation function
    function dydx = beamODE(x,y)
        dydx = [y(2); y(3); y(4); 0]; % Since EI*w'''' = 0
    end

    % Define the boundary conditions function to match analytical scenario
    function res = beamBC(ya,yb)
        % Adjust to reflect moment M0 applied at the free end
        % Clamped end: w(0) = 0, w'(0) = 0
        % Free end: EI*w''(L) = M0, EI*w'''(L) = 0
        res = [ya(1); ya(2); EI*yb(3) - M0; yb(4)];
    end

    % Initial mesh and initial guess for the solution
    xmesh = linspace(0, L, 100);
    solinit = bvpinit(xmesh, [0 0 0 0]);

    % Solve the boundary value problem
    sol = bvp4c(@beamODE, @beamBC, solinit);

    % Plot the solution
    x = linspace(0, L, 100);
    y = deval(sol, x);
    figure;
    plot(x, y(1,:), 'LineWidth', 2);
    xlabel('Position along the beam (m)');
    ylabel('Deflection (m)');
    title('Deflection of a Clamped-Free Beam with Applied Moment at Free End');
    grid on;
end