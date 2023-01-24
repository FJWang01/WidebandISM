

close all
clear
clc

%%
n = 10:15;
m = 10*ones(size(n));
t = -1:1e-3:1;

Pn = assoc_Legendre_poly(n, m, t);


figure;
plot(t, squeeze(Pn));
