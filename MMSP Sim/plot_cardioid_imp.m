
close all
clear
clc

%% Load data
load('test_cardioid_imp.mat');

%% Plot impulse responses at mic No.14
mic_idx = 17;

h_mics_single = h_mics(mic_idx, :);
h_mics_TDW_single = real(h_mics_TDW(mic_idx, :));

figure;
tiledlayout(3, 1);
nexttile;
plot(time_vec(1:numel(h_mics_single)), h_mics_single./max(h_mics_single));
hold on;
plot(time_vec, h_mics_TDW_single./max(h_mics_TDW_single));
xlim([0, 0.04]);
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('SMIR', 'TDW-SMIR');
grid on;

nexttile;
plot(time_vec(1:numel(h_mics_single)), h_mics_single./max(h_mics_single));
hold on;
plot(time_vec, h_mics_TDW_single./max(h_mics_TDW_single));
xlim([3.6e-3, 5.2e-3]);
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('SMIR', 'TDW-SMIR');
grid on;

nexttile;
plot(time_vec(1:numel(h_mics_single)), h_mics_single./max(h_mics_single));
hold on;
plot(time_vec, h_mics_TDW_single./max(h_mics_TDW_single));
xlim([6e-3, 19e-3]);
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('SMIR', 'TDW-SMIR');
grid on;

%% Plot all impulse responses
h_mics_TDW_real = real(h_mics_TDW);
%./max(h_mics, [], 'all')
f1 = figure;
tiledlayout(1, 2, 'TileSpacing', 'compact');
nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics./max(h_mics, [], 'all'));
xlim([0, 0.02]);
xlabel('Time (seconds)');
ylabel('Microphone index');
caxis([-0.1, 0.4]);

nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics_TDW_real(:, 1:size(h_mics, 2))./max(h_mics_TDW_real, [], 'all'));
xlim([0, 0.02]);
xlabel('Time (seconds)');
ylabel('Microphone index');
caxis([-0.1, 0.4]);
cB6 = colorbar;
cB6.Label.String = 'Amplitude';
cB6.Layout.Tile = 'east';

%% 
f1.PaperType = 'a4';
f1.PaperOrientation = 'landscape';
print('cardioid_imps_all_mics', '-dpdf', '-bestfit');