% This file plots the MTP filters for multiple degrees in n

close all
clear
clc

%% Set-up
[n, ~, ~] = getDegreeOrderPairs(5); % degrees of Legendre polynomial

% n = (0:5).';
fs = 96e3; % sampling frequency, in Hz
Ts = 1/fs; % sampling period, in seconds
r1 = 0.34; % small radius, in m
r2 = 2.2; % big radius, in m
c = 343; % speed of sound, in m/s
Tmax = 2*(r1+r2)/c; % maximum time of interest
t = 0:Ts:Tmax; % time sample vector

%% Calculate MTP filter of degree n
x = MTP_filter_all_degrees(n, r1, r2, c, t);

%% Plotting
figure;
plot(t, x);
xlabel('Time (s)');

