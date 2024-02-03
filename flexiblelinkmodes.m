% Beam properties
L = 1; % Length of the beam in meters
a = 6.35e-3; % Side length of the square cross-section in meters
I = a^4 / 12; % Moment of inertia
A = a^2; % Cross-sectional area in m^2
rho = 2700; % Density of aluminum in kg/m^3
E = 69e9; % Young's modulus in Pa

% Calculate EI and rho*A
EI = E * I;
rhoA = rho * A;

% Define function for characteristic equation (clamped-free beam)
charEq = @(k) cos(k*L).*cosh(k*L) + 1;

% Initialize variables to store modal frequencies
k = zeros(1,5);
omega = zeros(1,5);
f = zeros(1,5);

% Find the first five roots of the characteristic equation numerically
options = optimset('Display', 'off'); % Suppress fsolve output
for n = 1:5
    % Initial guess for nth root, improving guess range for each mode
    kGuess = (2*n-1)*pi/(2*L);
    k(n) = fsolve(@(k) charEq(k), kGuess, options);
    
    % Calculate omega (rad/s) and f (Hz) for each mode
    omega(n) = sqrt(EI/rhoA) * k(n)^2;
    f(n) = omega(n) / (2*pi);
end

% Display results
disp('Modal Frequencies:');
for n = 1:5
    fprintf('Mode %d: k = %.4f rad/m, \omega = %.4f rad/s, f = %.4f Hz\n', n, k(n), omega(n), f(n));
end

% Given data from the previous MATLAB code example (approximated for demonstration)
k = [1.87510407, 4.69409113, 7.85475744, 10.99554073, 14.13716839]; % Approximated roots k for each mode
f = [5.18, 32.5, 90.99, 178.31, 294.76]; % Frequencies in Hz for each mode

% Plotting
figure;
plot(k, f, 'o-', 'Color', 'blue', 'MarkerSize', 8, 'LineWidth', 2);
title('Modal Frequencies vs. Wave Numbers (k)');
xlabel('Wave Number k (rad/m)');
ylabel('Frequency f (Hz)');
grid on; % MATLAB automatically applies grid lines for both axes
