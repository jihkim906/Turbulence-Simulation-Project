clear; close all; clc;

%% Simulation Parameters
N = 2^8;                    % Grid points (power of 2 for FFT)
a = 0; b = 2*pi;           % Domain: [0, 2pi] × [0, 2pi]
L = abs(b-a);              % Domain length
h = L / N;                 % Grid spacing
step = h * (1:1:N);        % Grid coordinates
t0 = 0;                    % Initial time
T = 15;                    % Final time
viscosity = 0.001;         % Kinematic viscosity

%% Wavenumber Setup
wave_num = 1i * 2 * pi / L * [0:N/2-1 ,-N/2:-1];
[x, y] = meshgrid(step, step);
[kx, ky] = meshgrid(wave_num, wave_num);

%% Laplacian Operator
Lap_hat = kx.^2 + ky.^2;
k2 = Lap_hat; 
k2(1, 1) = 1;              % Avoid division by zero

%% Initialize Arrays
Re_prob = []; Re_max = []; Re_mean = [];

%% Video Setup
FrameRate = 10;
k = 0;
j = 1;

%% Initial Vorticity Field
% Create multiple Gaussian vortices
vorticity = zeros(N, N);
num_gaussians = 15;
spacing = L / (num_gaussians + 1);
sigma = 0.5;
rng(25);

for row = 1:num_gaussians
    for col = 1:num_gaussians
        xc = col * spacing;
        yc = row * spacing;
        sign = randi([0, 1]) * 2 - 1;     % Random sign
        vorticity = vorticity + sign * exp(-((x - xc).^2 + (y - yc).^2) / sigma^2);
    end
end

vort_hat = fft2(vorticity);

%% Video File
vidObj = VideoWriter('Simulation.avi');
open(vidObj);

%% Main Time Loop
while t0 < T
    
    % Compute stream function
    psi_hat = -vort_hat ./ k2;
    
    % Recover velocity field
    u = real(ifft2(ky .* psi_hat));
    v = real(ifft2(-kx .* psi_hat));
    
    % Compute vorticity gradients
    wx = real(ifft2(kx .* vort_hat));
    wy = real(ifft2(ky .* vort_hat));
    
    % Advection term
    VgradW = u .* wx + v .* wy;
    VgradW_hat = fft2(VgradW);
    
    % Velocity magnitude for CFL condition
    u_sqrt = sqrt(u.^2 + v.^2);
    
    % Adaptive time step
    dt = h/norm(u_sqrt,inf)*5;
    
    % Crank-Nicolson time integration
    vort_hat_next = 1 ./ (1 / dt - 0.5 * viscosity * Lap_hat) .* ...
                   ((1 / dt + 0.5 * viscosity * Lap_hat) .* vort_hat - VgradW_hat);
    
    % Visualization
    w = real(ifft2(vort_hat_next));
    imagesc(w); 
    shading interp; 
    colorbar; 
    colormap('jet');
    
    % Display time
    text(0.1, 0.9, ['Time: ', num2str(t0, '%.3f')], ...
         'Units', 'normalized', 'Color', 'k', 'FontSize', 18);
    
    % Record frame
    currFrame = getframe(gcf);
    writeVideo(vidObj, currFrame);
    
    % Update for next time step
    vort_hat = vort_hat_next;
    t0 = t0 + dt;
    
end

% Close video file
close(vidObj);