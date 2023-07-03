
close all
clear
clc

%% Load loudspeaker DIRs
load('Genelec_8020.mat'); % DIRs
load('elevs_meas.mat'); % elevations of measurement points
load('azims_meas.mat'); % azimuths of measurement points
load('r_meas.mat'); % radii of measurement points
load('DirPat_sampling_rate.mat'); % DIR sampling frequency

% Truncate DIR in time domain
DIR = DIR(:, 1:512);

%% Calculate the time-domain SH coefficients of the loudspeaker DIR
N_src = 5;
[n_src, m_src, ~] = getDegreeOrderPairs(N_src);
hnm_src = forwardSHT_mode_match(DIR, elevs_meas, azims_meas, N_src);
DIR_recon = inverseSHT(hnm_src, elevs_meas, azims_meas, n_src, m_src);

figure;
mesh(DIR);
title('Original loudspeaker DIR');

figure;
mesh(real(DIR_recon));
title('Reconstructed loudspeaker DIR')

%% Constants
% Sampling frequencies
meas_fs = fs; % for DIR measurements
TDW_SMIR_fs = fs; % for TDW-SMIR simulation

% Speed of sound
c = 343; % metres per second

%% Room setup
% Reflection order
Reflection_order = 2;

% Room Properties
Room.size = [4, 6, 3]; % (Lx, Ly, Lz) in metres
Room.b = [0.45, 0.7, 0.8, 0.5, 0.6, 0.75]; % wall reflection coefficients

%% Source setup
Source.GL_Os = [1, 3.5, 2.1]; % location of the source w.r.t the origin of the room 

Source.GL_SH_coeffs = hnm_src;

Source.tau = (0:size(DIR, 2)-1)/meas_fs;

%% Receiver setup
Receiver.GL_Or = [2.5, 3.5, 2.1]; % centre of the observation sphere w.r.t the origin of the room
Receiver.R = 0.042; % radius of the observation sphere, 0.042 for Eigenmike

% Microphone locations [azimuth, elevation] w.r.t. the centre of the Eigenmike
hom = getEigenmike;
azim_mics = hom.az;
elev_mics = hom.el;

mic_locs_on_Eigenmike = [azim_mics, elev_mics]; 

% SH truncation of the observed signal
f_max = 4000; % max frequency of interest, in Hz
k_max = 2*pi*f_max/c; % max wavenumber of interest
N_obv = ceil(exp(1)*k_max*Receiver.R/2);

[n_obv, m_obv, ~] = getDegreeOrderPairs(N_obv);

%% Other setups
% Low-pass filter
LPF_order = 256;
LPF_fn = f_max;

t_max = 0.1; % maximum time of interest, in seconds
time_vec = 0:1/TDW_SMIR_fs:t_max; % time sample vector, in seconds

%% Visualise setup 
% Cartesian coordinates of microphones on the Eigenmike expressed w.r.t. 
% the centre of the Eigenmike
[x_mics_offset, y_mics_offset, z_mics_offset] = mysph2cart(azim_mics, ...
    elev_mics, Receiver.R);

% Cartesian coordinates of microphones on the Eigenmike expressed w.r.t.
% the origin of the room
x_mics = Receiver.GL_Or(1) + x_mics_offset; % x coordinate
y_mics = Receiver.GL_Or(2) + y_mics_offset; % y coordinate
z_mics = Receiver.GL_Or(3) + z_mics_offset; % z coordinate

% microphone indexes
mic_indexes = (1:numel(azim_mics)).';
mic_indexes_str = num2str(mic_indexes); 
mic_indexes_cell = cellstr(mic_indexes_str);

% Visualise setup
figure;
scatter3(x_mics, y_mics, z_mics, 'ro'); % microphones on the Eigenmike
text(x_mics, y_mics, z_mics, mic_indexes_cell); % add microphone indexes
hold on
scatter3(Receiver.GL_Or(1), Receiver.GL_Or(2), Receiver.GL_Or(3), 'g^'); % Centre of the Eigenmike
hold on
scatter3(Source.GL_Os(1), Source.GL_Os(2), Source.GL_Os(3), 'b+'); % Centre of the source
xlim([0, Room.size(1)]);
ylim([0, Room.size(2)]);
zlim([0, Room.size(3)]);
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');

%% TDW-SMIR
Image = getImageSourcesISM_forLoop_SC(Reflection_order, Room, Source, Receiver, TDW_SMIR_fs, N_obv, c, time_vec, LPF_order, LPF_fn);

%% Calculate the observed signal using region-to-rgion RTF
RoomCC.b = Room.b;
RoomCC.size = Room.size;
RoomCC.D = Reflection_order;

SourceCC.Os = Source.GL_Os;
SourceCC.Rs = 1;
SourceCC.N = N_src;

ReceiverCC.Or = Receiver.GL_Or;
ReceiverCC.Rr = Receiver.R;
ReceiverCC.V = N_obv;

% f_vec = 0:f_max/128:f_max;

%% DIR frequency response
N_FFT = 4096;
DTF_full_band = fft(DIR, N_FFT, 2);
f_vec_full_band = 0:fs/N_FFT:fs;
f_vec = f_vec_full_band(f_vec_full_band <= f_max);
DTF = DTF_full_band(:, f_vec_full_band <= f_max);

% DTF = zeros(numel(elevs_meas), 129);
% for meas_idx = 1:numel(elevs_meas)
%     DTF(meas_idx, :) = freqz(DIR(meas_idx, :), 1, f_vec, fs);
% end

k_vec = 2*pi*f_vec/c;

%% Mode coupling matrix
A = zeros((ReceiverCC.V+1)^2, (SourceCC.N+1)^2, numel(k_vec));
for k_idx = 1:numel(k_vec)
    [Image_CC, A(:, :, k_idx)] = couplingcoefficients(k_vec(k_idx), RoomCC, SourceCC, ReceiverCC);
end
%% Obverved signal
% Source SH coefficients
Hnm_src_DTF = forwardSHT_mode_match(DTF, elevs_meas, azims_meas, SourceCC.N);

% Normalise using hn(krs)
hn_src_mat = sph_Hankel_1_adapted(n_src, k_vec, SourceCC.Rs);
Hnm_src_DTF_norm = Hnm_src_DTF./hn_src_mat;

% Observed signal's SH coefficients
Hnm_obv = zeros(size(A, 1), numel(k_vec));
for k_idx = 1:numel(k_vec)
    Hnm_obv(:, k_idx) = A(:, :, k_idx) * Hnm_src_DTF_norm(:, k_idx);
end

% Multiply by jn(kr)
jn_obv_mat = sph_Bessel_1_adapted(n_obv, k_vec, ReceiverCC.Rr);
Hnm_obv_SMIR = Hnm_obv.*jn_obv_mat;

%% Observed signal
h_mics_TDW_SMIR = inverseSHT(Image.Obv_GL_SH_coeff_sum, elev_mics, azim_mics, n_obv, m_obv);
figure;
mesh(real(h_mics_TDW_SMIR));

% H_mics_TDW_SMIR = zeros(numel(elev_mics), 129);
% for mic_idx = 1:numel(elev_mics)
%     H_mics_TDW_SMIR(mic_idx, :) = freqz(h_mics_TDW_SMIR(mic_idx, :), 1, f_vec, fs ); 
% end

H_mics_TDW_SMIR_full_band = fft(h_mics_TDW_SMIR, N_FFT, 2);
H_mics_TDW_SMIR = H_mics_TDW_SMIR_full_band(:, f_vec_full_band <= f_max);

H_mics_SMIR = inverseSHT(Hnm_obv_SMIR, elev_mics, azim_mics, n_obv, m_obv);


%% Plot TDW-SMIR results in the time domain
% for mic_idx = 1:numel(elev_mics)
%     figure;
%     plot(real(h_mics_TDW_SMIR(mic_idx, :)));
%     % saveas(gcf,['Signal at mic No.',  num2str(mic_idx),'.fig']);
% end

%% Compare results in the frequency domain
H_mics_TDW_SMIR_max = max(mag2db(abs(H_mics_TDW_SMIR)), [], 'all');
H_mics_SMIR_max = max(mag2db(abs(H_mics_SMIR)), [], 'all');

H_mics_TDW_SMIR_ref = mean(mag2db(abs(H_mics_TDW_SMIR)), 'all', 'omitnan');
H_mics_SMIR_ref = mean(mag2db(abs(H_mics_SMIR)), 'all', 'omitnan');


for mic_idx = 1:numel(elev_mics)
    H_mics_TDW_SMIR_single = H_mics_TDW_SMIR(mic_idx, :);
    H_mics_SMIR_single = H_mics_SMIR(mic_idx, :);
    figure;
    plot(f_vec, mag2db(abs(H_mics_TDW_SMIR_single)) - H_mics_TDW_SMIR_max);
    hold on;
    plot(f_vec, mag2db(abs(H_mics_SMIR_single)) - H_mics_SMIR_max);
    xlim([100, 4e3]);
    xticks([100, 500:500:4e3]);
end

