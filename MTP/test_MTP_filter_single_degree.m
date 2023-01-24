% This file plots the MTP filter for a single degree n

close all
clear
clc

%% Set-up
n = 10; % degree of Legendre polynomial
fs = 44100; % sampling frequency, in Hz
Ts = 1/fs; % sampling period, in seconds
r1 = 0.5; % small radius, in m
r2 = 2; % big radius, in m
c = 343; % speed of sound, in m/s
Tmax = 2*(r1+r2)/c; % maximum time of interest
t = 0:Ts:Tmax; % time sample vector

%% Calculate MTP filter of degree n
x = MTP_filter_single_degree(n, r1, r2, c, t);

%% Plotting
figure;
plot(t, x);
xlabel('Time (s)');




