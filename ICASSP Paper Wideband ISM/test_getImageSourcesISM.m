
close all
clear
clc

%% Set-up
% Sampling frequency
fs = 16e3; % Hz

% Speed of sound
c = 343; % metres per second

% Reflection order
Reflection_order = 2;

% Room Properties
Room.size = [4, 6, 3]; % (Lx, Ly, Lz) in metres
Room.b = [0.45, 0.7, 0.8, 0.5, 0.6, 0.75]; % wall reflection coefficients

% Source
Source.GL_Os = [1.5, 3.4, 2.4]; % location of the source w.r.t the origin of the room 
Source.GL_SH_coeffs = [1, 0;0, 0;0, 1;0,0]; % degree of source directivity pattern
% Source.GL_SH_coeffs = rand(16,2);
Source.tau = [0, 0.01];

% Receiver
Receiver.GL_Or = [1.5, 3.4, 1]; % location of the receiver w.r.t the origin of the room
Receiver.R = 0.2; % radius of the observation sphere

% Spherical harmonic truncation order
f_max = 2000; % max frequency of interest, in Hz
k_max = 2*pi*f_max/c; % max wavenumber of interest
N = ceil(exp(1)*k_max*Receiver.R/2);

% Low-pass filter
LPF_order = 256;
LPF_fn = f_max;

t_max = 0.2; % maximum time of interest, in seconds
time_vec = 0:1/fs:t_max; % time sample vector, in seconds
%% Find all image sources
Image = getImageSourcesISM(Reflection_order, Room, Source, Receiver, fs, N, c, time_vec, LPF_order, LPF_fn);

%% Reconstruction grid
[n, m, ~] = getDegreeOrderPairs(N);

elev_samples = 0:pi/50:pi; % elevation samples
azim_samples = 0:pi/50:2*pi; % azimuth samples at each elevation sample 

[elev_grid, azim_grid] = meshgrid(elev_samples, azim_samples); % construct 2D grid

elev_grid_col = elev_grid(:);
azim_grid_col = azim_grid(:);

p_grid_col = inverseSHT(Image.Obv_GL_SH_coeff_final, elev_grid_col, azim_grid_col, n, m);

% Reshape p_grid_col
p_grid = zeros(size(elev_grid, 1), size(elev_grid, 2), numel(time_vec));
for t_idx = 1:numel(time_vec)
    p_grid(:, :, t_idx) = reshape(p_grid_col(:, t_idx), size(elev_grid));
end
%% Plot as a function of azimuth
elev_idx = 26;
elev_plot = rad2deg(elev_samples(elev_idx)); % chosen elevation angle in degree

% size(p_grid_at_elev) = [numel(azim_samples), numel(time_vec)]
p_grid_at_elev = squeeze(p_grid(:, elev_idx, :)); 

figure;
imagesc(time_vec, azim_samples, real(p_grid_at_elev));
xlabel('Time (s)');
ylabel('Azimuth angle (radians)');
zlabel('Real part of amplitude');
title(['Real part of the reconstructed signal at elevation ', num2str(elev_plot), ' degrees']);
xlim([0, 0.04]);
caxis([-80, 60]);

