
% This file compares the open SMA signals due to a source whose frequency
% invariant directivity is a cardioid
close all
clear
clc

%% Constants
% Sampling frequency
fs = 44.1e3; % Hz

% Speed of sound
c = 343; % metres per second

% Density of air at 20 degrees Celcius
rho = 1.2041; % density of air, in kg/m^3 at 20 degrees Celcius

%% Room setup
% Reflection order
Reflection_order = 2;

% Room Properties
Room.size = [4, 6, 3]; % (Lx, Ly, Lz) in metres
Room.b = [0.45, 0.7, 0.8, 0.5, 0.6, 0.75]; % wall reflection coefficients

% Room property specific to SMIR generator
refl_coeff_ang_dep = 0; % Real reflection coeff(0) or angle dependent reflection coeff(1)

%% Source setup 
Source.GL_Os = [1, 3.5, 2.1]; % location of the source w.r.t the origin of the room

Source.tau = 0;

% Source setup specific to SMIR generator
% Cardioid
src_type = 'c'; % Source directivity pattern ('o','c','s','h','b')

%% Receiver setup
Receiver.GL_Or = [2.5, 3.5, 2.1]; % centre of the observation sphere w.r.t the origin of the room
Receiver.R = 0.042; % radius of the observation sphere, 0.042 for Eigenmike

% The following three lines are for higher-order source
[src_ang(1),src_ang(2)] = mycart2sph(Receiver.GL_Or(1)-Source.GL_Os(1),...
    Receiver.GL_Or(2)-Source.GL_Os(2),Receiver.GL_Or(3)-Source.GL_Os(3)); 
SH_rot_mat = WignerDM(1, src_ang(1), src_ang(2), 0);

% Source SH coefficients for a cardioid
Source.GL_SH_coeffs = SH_rot_mat * [sqrt(4*pi);0;sqrt(4*pi/3);0];

% Receiver setup specific to SMIR Generator
sphType = 'open';  % Type of sphere (open/rigid), keep 'open'

%% Microphone locations [azimuth, elevation] w.r.t. the origin of the Eigenmike
hom = getEigenmike;
azim_mics = hom.az;
elev_mics = hom.el;

mic_locs_on_Eigenmike = [azim_mics, elev_mics]; 

%% Spherical harmonic truncation order
f_max = 4000; % max frequency of interest, in Hz
k_max = 2*pi*f_max/c; % max wavenumber of interest
N = ceil(exp(1)*k_max*Receiver.R/2);

%% Other setups
num_IR_samples = 2*1024; % number of samples of the impulse response
OversamplingFactor = 2; % Oversampling factor

HP = 0; % Enable 4-th order Butterworth high pass filter with cut-off frequency 50 Hz

% Low-pass filter
LPF_order = 512;
LPF_fn = f_max;

t_max = 0.2; % maximum time of interest, in seconds
time_vec = 0:1/fs:t_max; % time sample vector, in seconds

%% Plot setup
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

%% TDW - SMIR
Image = getImageSourcesISM(Reflection_order, Room, Source, Receiver, fs, N, c, time_vec, LPF_order, LPF_fn);

%% Original SMIR
[h_mics, H_mics] = smir_generator(c, fs, Receiver.GL_Or, Source.GL_Os, Room.size, Room.b, sphType, ...
    Receiver.R, mic_locs_on_Eigenmike, N, num_IR_samples, OversamplingFactor, Reflection_order, ...
    refl_coeff_ang_dep, HP, src_type, src_ang);

%% Plot observed signals
[n, m, ~] = getDegreeOrderPairs(N);
h_mics_TDW = inverseSHT(Image.Obv_GL_SH_coeff_final, elev_mics, azim_mics, n, m);

figure;
mesh(real(h_mics_TDW));

max_h_mics_TDW = max(real(h_mics_TDW), [], 'all');

figure;
mesh(h_mics);

max_h_mics = max(h_mics, [], 'all');


%% Comparison 
for mic_idx = 1:numel(elev_mics)
    figure;
    plot(real(h_mics_TDW(mic_idx, :))./max_h_mics_TDW);
    hold on;
    plot(h_mics(mic_idx, :)./max_h_mics);
    xlim([1, 2049]);
    legend('ICASSP', 'SMIR');
    %saveas(gcf,['Signal at mic No.',  num2str(mic_idx),'.fig']);
end


