% Beam properties
L = 1;                   % Length of the beam in meters
a = 6.35e-3;             % Side length of the square cross-section in meters
I = a^4 / 12;            % Moment of inertia in m^4
A = a^2;                 % Cross-sectional area in m^2
rho = 2700;              % Density of aluminum in kg/m^3
E = 69e9;                % Young's modulus in Pa

% Calculate EI and rho*A
EI = E * I;
rhoA = rho * A;

% Parameters
nx = 15;                 % Length of space domain
nt = 10^5;               % Length of time domain
tf = 1;                 % Time of simulation in seconds

% Compute the mesh spacing and time step
dx = L / nx;
dt = tf / nt;

% Drawing coordinates
x = linspace(0, L, nx + 1); 
t = linspace(0, tf, nt + 1);

% Create memory to save data
w = zeros(nx + 1, nt + 1); 

% Initial Conditions
P = -5 * 9.81; % Tip payload in Newtons
for i = 1:nx + 1
    w(i, 1) = -(P * L / (EI)) * ((x(i)^2) / 2 - (x(i)^3) / (6 * L));
    w(i, 2) = w(i, 1);
end

% Damping, payload, and other parameters
ksi = 0.2;    % Damping ratio
Ms = 5;       % Mass of the tip payload in kg
alpha = (dt * 2 * ksi * sqrt(EI / rhoA) / (rhoA));

% Boundary disturbance
d = 0.1 + 0.1 * sin(pi * (0:nt) * dt) + 0.1 * sin(2 * pi * (0:nt) * dt) + 0.1 * sin(3 * pi * (0:nt) * dt);

% Main loop for FDM
for j = 3:nt + 1     
    for i = 3:nx - 1
        wxxxx = (w(i + 2, j - 1) - 4 * w(i + 1, j - 1) + 6 * w(i, j - 1) - 4 * w(i - 1, j - 1) + w(i - 2, j - 1)) / dx^4;
        w(i, j) = ((2 + alpha) * w(i, j - 1) - w(i, j - 2) + (-EI * wxxxx) * (dt^2 / rhoA)) / (1 + alpha);
    end

    % Boundary Conditions
    wxL = (w(nx + 1, j - 1) - w(nx, j - 1)) / dx;
    wxxxL = (w(nx + 1, j - 1) - 3 * w(nx, j - 1) + 3 * w(nx - 1, j - 1) - w(nx - 2, j - 1)) / dx^3;
    w(nx + 1, j) = 2 * w(nx + 1, j - 1) - w(nx + 1, j - 2) + (d(j - 1) + EI * wxxxL) * dt^2 / Ms;
    w(nx, j) = (w(nx + 1, j) + w(nx - 1, j)) / 2; 
end

% Visualization (Mesh plot)
figure(1)
mesh(x, t, w');
title('Displacement of Beam with Boundary and Initial Conditions');
ylabel('Time (s)', 'FontSize', 12);
xlabel('x (m)', 'FontSize', 12);
zlabel('Displacement (m)', 'FontSize', 12);
view([60 45]);

% Open a video file
videoFileName = 'flexiblelinkvibration.avi';
v = VideoWriter(videoFileName, 'Motion JPEG AVI'); % Specify the format for broader compatibility
v.FrameRate = 60; % Increase the frame rate for a smoother and faster playback
open(v);

% Visualization (Faster Animation with recording)
figure;
skipFrames = 100; % Define how many frames to skip for a faster animation
for j = 1:skipFrames:length(t) % Update loop to increment by skipFrames
    plot(x, w(:, j));
    xlabel('x (m)');
    ylabel('Displacement (m)');
    title(sprintf('Beam Displacement at Time = %.2f s', t(j)));
    axis([0 L min(w(:)) max(w(:))]);
    drawnow;
    
    % Capture the plot as an image 
    frame = getframe(gcf);
    writeVideo(v, frame); % Write frame to video
    
    % Removed pause for a faster execution
end

% Close the video file
close(v);

% Notify user
fprintf('Fast animation has been recorded and saved as %s\n', videoFileName);