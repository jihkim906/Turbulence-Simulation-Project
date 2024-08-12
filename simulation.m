clear; close all; clc;

% Simulation Setting
N = 2^8; a = 0; b = 2*pi; L = abs(b-a);
h = L / N; step = h * (1:1:N); t0 = 0; T = 15; viscosity = 0.001;

% Wave Number
wave_num = 1i * 2 * pi / L * [0:N/2-1 ,-N/2:-1];

% Grid Setting
[x, y] = meshgrid(step, step);
[kx, ky] = meshgrid(wave_num, wave_num);

Lap_hat = kx.^2 + ky.^2;
k2 = Lap_hat; k2(1, 1) = 1;

% Initialize Reynolds Number Lists
Re_prob = []; Re_max = []; Re_mean = [];


% Movie File Data Allocation Setup
FrameRate = 10;
k = 0;
j = 1;

% Initial Vorticity Setting
vorticity = zeros(N, N);

% Parameters for Gaussian functions
num_gaussians = 15;
spacing = L / (num_gaussians + 1);
sigma = 0.5; % Standard deviation for the Gaussians

% Set random seed for reproducibility
rng(25);

% Loop to add Gaussians
for row = 1:num_gaussians
    for col = 1:num_gaussians
        xc = col * spacing;
        yc = row * spacing;
        sign = randi([0, 1]) * 2 - 1; % Random sign: +1 or -1
        
        vorticity = vorticity + sign * exp(-((x - xc).^2 + (y - yc).^2) / sigma^2);
    end
end

% FFT transformation of initial vorticity field
vort_hat = fft2(vorticity);

% Prepare the new file.
vidObj = VideoWriter('Simulation.avi');
open(vidObj);

% Method
while t0 < T
    
    % stream function and velocity
    psi_hat = -vort_hat ./ k2;
    u = real(ifft2(ky .* psi_hat));
    v = real(ifft2(-kx .* psi_hat));
    wx = real(ifft2(kx .* vort_hat));
    wy = real(ifft2(ky .* vort_hat));
    
    % VgradW calculation
    VgradW = u .* wx + v .* wy;
    VgradW_hat = fft2(VgradW);

    % reynolds number calculation
    u_sqrt = sqrt(u.^2+v.^2);
    rey_prob = norm(VgradW)/norm(viscosity*real(ifft2(k2 .* vort_hat)));
    Re_prob = [Re_prob, rey_prob];
    rey_max = norm(u_sqrt,Inf)*sqrt(sigma) / viscosity;
    Re_max = [Re_max, rey_max];
    rey_mean = mean(mean(u_sqrt))*sqrt(sigma) / viscosity;
    Re_mean = [Re_mean, rey_mean];

    % adaptive time step
    dt = h/norm(u_sqrt,inf)*5;
 
    % Crank-Nicholson
    vort_hat_next = 1 ./ (1 / dt - 0.5 * viscosity * Lap_hat) .* ((1 / dt + 0.5 * viscosity * Lap_hat) .* vort_hat - VgradW_hat);

    % plot vorticity field
    w = real(ifft2(vort_hat_next));
    imagesc(w); shading interp; colorbar; colormap('jet'); 
    % Display time
    text(0.1, 0.9, ['Time: ', num2str(t0, '%.3f')], 'Units', 'normalized', 'Color', 'k', 'FontSize', 18);

     % Write each frame to the file.
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);

    % update vorticity field
    vort_hat = vort_hat_next;

    % update time
    t0 = t0 + dt;

   
end

% Close the file
close(vidObj);

% Plot Re1, Re2, Re3 in separate subplots
figure;
subplot(3, 1, 1);
plot(Re_prob, 'r');
xlabel('Time Step');
ylabel('Re1');
title('Reynolds number based on the vorticity gradient');

subplot(3, 1, 2);
semilogy(Re_max, 'g');
xlabel('Time Step');
ylabel('Re2');
title('Maximum Reynolds number');

subplot(3, 1, 3);
semilogy(Re_mean, 'b');
xlabel('Time Step');
ylabel('Re3');
title('Mean Reynolds number');
