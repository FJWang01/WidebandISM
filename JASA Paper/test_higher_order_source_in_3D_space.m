
close all
clear
clc

%% Set-up
r = 0.2; % radius of the observation sphere, in metres
r_s = 1.5; % radius of the source, in metres
elev_s = pi/2; % elevation angle of the source, in radian
azim_s = 2*pi/3; % azimuth angle of the source, in radian

c = 343; % speed of sound, in metres per second
fs = 48e3; % sampling frequency, in Hz

up_factor = 3; % unsampling factor
fup = fs * up_factor;  % sampling frequency after upsampling

f_max = 2.5e3; % max frequency of interest, in Hz
k_max = 2*pi*f_max/c; % max wavenumber corresponding to max frequency of interest
N = ceil(exp(1)*k_max*r/2); % spherical harmonic truncation order

v = 1; % degree of the source pattern
u = 0; % order of the source pattern

%% Find the time sample vector
t_max = 4*(r+r_s)/c; % max time of interest, in seconds
time_vec_up = 0:1/fup:t_max; % time sample vector

time_vec = downsample(time_vec_up, up_factor);

%% Find SH coefficients
% size(SH_coeff_up) = [(N+1)^2, numel(time_vec_up)]
SH_coeff_up = higher_order_source_in_3D_space(r, r_s, elev_s, azim_s, v, u, N, time_vec_up, c);

%% Design LPF to extract information between [0, fs/2]
[b, a] = cheby2(16, 60, fs/fup);

% b = fir1(255, fs/fup);

% Zero-phase filtering
SH_coeff_lp_up = filtfilt(b, a, SH_coeff_up.');

%% Downsample SH_coeff_lp_up so that the sampling frequency is back to fs
SH_coeff_lp_fs = downsample(SH_coeff_lp_up, up_factor);

SH_coeff_lp_fs = SH_coeff_lp_fs.';

%% Reconstruction grid
[n, m, ~] = getDegreeOrderPairs(N);

elev_samples = 0:0.05:pi; % elevation samples
azim_samples = 0:0.05:2*pi; % azimuth samples 

[elev_grid, azim_grid] = meshgrid(elev_samples, azim_samples); % construct 2D grid

elev_grid_col = elev_grid(:);
azim_grid_col = azim_grid(:);

p_grid_col = inverseSHT(SH_coeff_lp_fs, elev_grid_col, azim_grid_col, n, m);

% Reshape p_grid_col
p_grid = zeros(size(elev_grid, 1), size(elev_grid, 2), numel(time_vec));
for t_idx = 1:numel(time_vec)
    p_grid(:, :, t_idx) = reshape(p_grid_col(:, t_idx), size(elev_grid));
end
%% Plot as a function of azimuth
elev_idx = 32;
elev_plot = rad2deg(elev_samples(elev_idx)); % chosen elevation angle in degree

% size(p_grid_at_elev) = [numel(azim_samples), numel(time_vec)]
p_grid_at_elev = squeeze(p_grid(:, elev_idx, :)); 

figure;
mesh(time_vec, azim_samples, real(p_grid_at_elev));
xlabel('Time (s)');
ylabel('Azimuth angle (radians)');
zlabel('Real part of amplitude');
title(['Real part of the reconstructed signal at elevation ', num2str(elev_plot), ' degrees']);

figure;
mesh(time_vec, azim_samples, imag(p_grid_at_elev));
xlabel('Time (s)');
ylabel('Azimuth angle (radians)');
zlabel('Imaginary part of amplitude');
title(['Imaginary part of the reconstructed signal at elevation ', num2str(elev_plot), ' degrees']);
