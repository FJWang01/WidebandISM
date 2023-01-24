% This file plots the MTP filter for plane wave

close all
clear
clc

%% Set-up
% [n, ~, ~] = getDegreeOrderPairs(5);

n = (0:5).'; % degree of Legendre polynomial
fs = 44100; % sampling frequency, in Hz
Ts = 1/fs; % sampling period, in seconds
r = 0.5; % small radius, in m
c = 343; % speed of sound, in m/s
Tmax = 2*r/c; % maximum time of interest
t = -Tmax:Ts:Tmax; % time sample vector

%% Calculate PW MTP filters of all degrees in n
x = PW_MTP_filter_all_degrees(n, r, c, t);

%% Plotting
figure;
plot(t, x);
xlabel('Time (s)');
ylabel('Amplitude');