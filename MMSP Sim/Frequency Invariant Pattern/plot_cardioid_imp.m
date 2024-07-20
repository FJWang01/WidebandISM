
close all
clear
clc

%% Load data
load('test_cardioid_imp.mat');

%% Plot impulse responses at mic No.17
mic_idx = 17;

h_mics_single = h_mics(mic_idx, :);
h_mics_TDW_single = real(h_mics_TDW(mic_idx, :));

figure;
tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','tight');
nexttile([2,1]);
plot(time_vec(1:numel(h_mics_single)), h_mics_single./max(h_mics_single));
hold on;
plot(time_vec, h_mics_TDW_single./max(h_mics_TDW_single));
xlim([0, 0.025]);
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
ylim([-0.05, 0.35]);
xlabel('Time (seconds)');
ylabel('Amplitude');
legend('SMIR', 'TDW-SMIR');
grid on;

%% Plot all impulse responses
% Define colormap
% mymap = cmap_redblue(256);
mymap1 = cmap_redwhiteblue(-0.05, 1, 256);

mymap2 = cmap_redwhiteblue(-0.05, 0.2, 256);

h_mics_TDW_real = real(h_mics_TDW);
%./max(h_mics, [], 'all')
f1 = figure;
tiledlayout(2, 2, 'TileSpacing', 'compact');

ax1 = nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics./max(h_mics, [], 'all'));
xlim([0, 0.01]);
xlabel('Time (seconds)');
ylabel('Microphone index');
colormap(ax1, mymap1);
clim([-0.05, 1]);

ax2 = nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics_TDW_real(:, 1:size(h_mics, 2))./max(h_mics_TDW_real, [], 'all'));
xlim([0, 0.01]);
xlabel('Time (seconds)');
ylabel('Microphone index');
clim([-0.05, 1]);
colormap(ax2, mymap1);
cB6 = colorbar;
cB6.Label.String = 'Amplitude';
cB6.Layout.Tile = 'east';

ax3 = nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics./max(h_mics, [], 'all'));
xlim([0.01, 0.02]);
xlabel('Time (seconds)');
ylabel('Microphone index');
colormap(ax3, mymap2);
clim([-0.05, 0.2]);

ax4 = nexttile;
imagesc(time_vec(1:size(h_mics, 2)), 1:32, h_mics_TDW_real(:, 1:size(h_mics, 2))./max(h_mics_TDW_real, [], 'all'));
xlim([0.01, 0.02]);
xlabel('Time (seconds)');
ylabel('Microphone index');
clim([-0.05, 0.2]);
colormap(ax4, mymap2);
cB7 = colorbar;
cB7.Label.String = 'Amplitude';
cB7.Layout.Tile = 'east';

%% Plot all impulse responses - waterfall plot
% Get label for each plot
for mic_id = 1:32
    mic_label{mic_id} = strcat('No.', num2str(mic_id), ...
        ' (', num2str(elev_mics(mic_id), 3), ', ',num2str(azim_mics(mic_id), 3), ') rad');
end

figure;
tiledlayout(8, 4, 'TileSpacing','tight', 'Padding','tight');
for mic_id = 1:32
    h_mics_single = h_mics(mic_id, :);
    h_mics_TDW_single = real(h_mics_TDW(mic_id, :));
    nexttile;
    plot(time_vec(1:numel(h_mics_single)), h_mics_single./max(h_mics, [], 'all'));
    hold on;
    plot(time_vec, h_mics_TDW_single./max(h_mics_TDW_real, [], 'all'));
    xlim([0, 0.025]);
    ylim([-0.1, 1]);
    text(0.005,0.8,mic_label(mic_id));
    if rem(mic_id, 4) == 1
        ylabel('Amplitude');
    end
    if (mic_id >= 29) && (mic_id <= 32)
        xlabel('Time (seconds)');
    end
    %legend('SMIR', 'TDW-SMIR');
    grid on;
end