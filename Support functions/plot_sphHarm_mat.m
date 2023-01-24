function Ynm = plot_sphHarm_mat(n, m, num_samples)
% This function plots the spherical harmonic function of 
% degree n and order m
%
% Input
% n - degree, scalar
% m - order, scalar
% num_samples - number of samples on the unit sphere, scalar
%
% Output
% Ynm - size(Ynm) = [1, 1, numel(num_samples)]

%% Generate sampling points on unit sphere
[X, Y, Z] = sphere(num_samples);
[azim, elev, ~] = cart2sph(X, Y, Z);
elev = pi/2 - elev; 

azim_col = azim(:);
elev_col = elev(:);
%% Convert to Cartesian coordinate and plot as desired
Ynm_3D = sphHarm_mat(n, m, elev_col, azim_col);
Ynm_col = squeeze(Ynm_3D);
Ynm = reshape(Ynm_col, size(elev));

%% Plot real part 
Ynm_real = real(Ynm);    
figure;
% plot positive values
Ynm_real_pos = zeros(size(Ynm_real));
Ynm_real_pos(Ynm_real>=0) = Ynm_real(Ynm_real>=0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_real_pos);
mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
hold on
% plot negative values
Ynm_real_neg = zeros(size(Ynm_real));
Ynm_real_neg(Ynm_real<0) = Ynm_real(Ynm_real<0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_real_neg));
mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
title(['Real part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
legend('Postive value', 'Negative value');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%% Plot imaginary part 
Ynm_imag = imag(Ynm);    
figure;
% plot positive values
Ynm_imag_pos = zeros(size(Ynm_imag));
Ynm_imag_pos(Ynm_imag>=0) = Ynm_imag(Ynm_imag>=0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_imag_pos);
mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
hold on
% plot negative values
Ynm_imag_neg = zeros(size(Ynm_imag));
Ynm_imag_neg(Ynm_imag<0) = Ynm_imag(Ynm_imag<0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_imag_neg));
mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
title(['Imaginary part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
legend('Postive value', 'Negative value');
grid on;
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

end