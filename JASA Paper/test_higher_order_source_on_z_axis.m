
close all
clear
clc

%% Set-up
r = 0.2; % radius of the observation sphere, in metres
r_s = 1.5; % radius of the source, in metres

c = 343; % speed of sound, in metres per second
fs = 48e3; % sampling frequency, in Hz

up_factor = 3; % upsampling factor
fup = fs * up_factor; % sampling frequency after upsampling

f_max = 20e3; % max frequency of interest, in Hz
k_max = 2*pi*f_max/c; % max wavenumber corresponding to max frequency of interest
% N = ceil(exp(1)*k_max*r/2); % spherical harmonic truncation order
N = 30;

v = 1; % degree of the source pattern
u = 0; % order of the source pattern

%% Find the time sample vector
t_max = 4*(r+r_s)/c; % max time of interest, in seconds
time_vec_up = 0:1/fup:t_max; % time sample vector for upsampled signal

time_vec = downsample(time_vec_up, up_factor);
%% Find SH coefficients
% size(SH_coeff_up) = [(N+1)^2, numel(time_vec_up)]
SH_coeff_up = higher_order_source_on_z_axis(r, r_s, v, u, N, time_vec_up, c);

%% Design LPF to extract information between [0, fs/2]
[b, a] = cheby2(8, 40, fs/fup);

% b = fir1(255, fs/fup);

% Zero-phase filtering
SH_coeff_lp_up = filtfilt(b, a, SH_coeff_up.');

%% Downsample SH_coeff_lp_up so that the sampling frequency is back to fs
SH_coeff_lp_fs = downsample(SH_coeff_lp_up, up_factor);

SH_coeff_lp_fs = SH_coeff_lp_fs.';
%% Reconstruction grid
[n, m, ~] = getDegreeOrderPairs(N);

elev_samples = 0:pi/100:pi; % elevation samples
azim_samples = 0:pi/100:2*pi; % azimuth samples 

[elev_grid, azim_grid] = meshgrid(elev_samples, azim_samples); % construct 2D grid

elev_grid_col = elev_grid(:);
azim_grid_col = azim_grid(:);

p_grid_col = inverseSHT(SH_coeff_lp_fs, elev_grid_col, azim_grid_col, n, m);

p_grid = zeros(size(elev_grid, 1), size(elev_grid, 2), numel(time_vec));
for t_idx = 1:numel(time_vec)
    p_grid(:, :, t_idx) = reshape(p_grid_col(:, t_idx), size(elev_grid));
end

%% Plot as a function of elevation
% Choose an azimuth angle
azim_plot_idx = 1; 
azim_plot = rad2deg(azim_samples(azim_plot_idx)); % chosen azimuth in degree

% size(p_recon_at_azim) = [numel(elev_samples), numel(time_vec)]
p_grid_at_azim = squeeze(p_grid(azim_plot_idx, :, :)); 

% Select three elevation indexes
elev_idx_1 = 26;
elev_idx_2 = 51;
elev_idx_3 = 76;

elev_at_idx_1 = elev_samples(elev_idx_1);
elev_at_idx_2 = elev_samples(elev_idx_2);
elev_at_idx_3 = elev_samples(elev_idx_3);

figure;
subplot(2, 2, 1);
imagesc(time_vec, elev_samples, real(p_grid_at_azim));
xlim([0, 6e-3]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Colatitude $\theta$ (radians)', 'interpreter', 'latex');
cBar = colorbar;
cBar.Label.String = 'Amplitude';
caxis([-1e3, 200]);
title('(a)');

subplot(2, 2, 2);
plot(time_vec, real(p_grid_at_azim(elev_idx_1, :)), ...
    'DisplayName',['\theta = ', num2str(elev_at_idx_1), ' rad']);
hold on
plot(time_vec, real(p_grid_at_azim(elev_idx_2, :)), ...
    'DisplayName',['\theta = ', num2str(elev_at_idx_2), ' rad']);
hold on
plot(time_vec, real(p_grid_at_azim(elev_idx_3, :)), ...
    'DisplayName',['\theta = ', num2str(elev_at_idx_3), ' rad']);
xlim([3.5e-3, 5.5e-3]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Amplitude', 'interpreter', 'latex');
legend;
title('(b)');


%% Plot as a function of azimuth
% Choose an elevation angle
elev_plot_idx = 51;
elev_plot = rad2deg(elev_samples(elev_plot_idx)); % chosen elevation in degree

% size(p_grid_at_elev) = [numel(azim_samples), numel(time_vec)]
p_grid_at_elev = squeeze(p_grid(:, elev_plot_idx, :)); 

azim_idx_1 = 51;
azim_idx_2 = 101;
azim_idx_3 = 151;

azim_at_idx_1 = azim_samples(azim_idx_1);
azim_at_idx_2 = azim_samples(azim_idx_2);
azim_at_idx_3 = azim_samples(azim_idx_3);


subplot(2, 2, 3);
imagesc(time_vec, azim_samples, real(p_grid_at_elev));
xlim([0, 6e-3]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Azimuth $\phi$ (radians)', 'interpreter', 'latex');
cBar = colorbar;
cBar.Label.String = 'Amplitude';
caxis([-400, 50]);
title('(c)');


subplot(2, 2, 4);
plot(time_vec, real(p_grid_at_elev(azim_idx_1, :)), '-o',...
    'DisplayName',['\phi = ', num2str(azim_at_idx_1), ' rad']);
hold on
% plot(time_vec, real(p_grid_at_elev(azim_idx_2, :)), '-x',...
%     'DisplayName',['\phi = ', num2str(azim_at_idx_2), ' rad']);
% hold on
plot(time_vec, real(p_grid_at_elev(azim_idx_3, :)), '-*',...
    'DisplayName',['\phi = ', num2str(azim_at_idx_3), ' rad']);
xlim([3.5e-3, 5.5e-3]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Amplitude', 'interpreter', 'latex');
legend;
title('(d)');

%% 3D plots
figure;
mesh(time_vec, elev_samples, real(p_grid_at_azim));
xlim([0, time_vec(round(numel(time_vec)/2))]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Elevation angle (radians)', 'interpreter', 'latex');
%zlabel('Real part of amplitude');
title(['Real part of the reconstructed signal at azimuth ', num2str(azim_plot), ' degrees']);
%cBar = colorbar;
%cBar.Label.String = 'Amplitude';

figure;
mesh(time_vec, elev_samples, imag(p_grid_at_azim));
xlim([0, time_vec(round(numel(time_vec)/2))]);
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Elevation angle (radians)', 'interpreter', 'latex');
%zlabel('Imaginary part of amplitude');
TT = title(['Imaginary part of the reconstructed signal at azimuth ', num2str(azim_plot), ' degrees']);
%cBar = colorbar;
%cBar.Label.String = 'Amplitude';