
function Ynm_rotated = plot_rotatedSphHarm(n, m, alpha, beta, gamma)
% This function plots the rotated spherical harmonic function of degree n
% and order m

% Input
% n - degree, scalar
% m - order, scalar
% alpha, beta, gamma - rotation angles

% Output
% Ynm_rotated - values of the rotated spherical harmonic function of degree
% n and order m

% Note:
% The rotation follows the zyz convention
% First, rotate gamma radian about z-axis
% Next, rotate beta radian about y-axia
% Then, rotate alpha radian about z-axis
%% Plotting grid
[X, Y, Z] = sphere(128);
[azim, elev, ~] = cart2sph(X, Y, Z);
elev = pi/2 - elev; 

%% Get rotated spherical harmonics
Ynm_rotated = getRotatedSphHarm(n, m, alpha, beta, gamma, elev, azim);

%% Plotting
%%%%%%%%%%%%%%%%%%%% Plot real part %%%%%%%%%%%%%%%%%%%%%%
Ynm_rotated_real = real(Ynm_rotated);    
figure;
% plot positive values
Ynm_real_pos = Ynm_rotated_real.*(Ynm_rotated_real>=0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_real_pos);
mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
hold on
% plot negative values
Ynm_real_neg = Ynm_rotated_real.*(Ynm_rotated_real<0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_real_neg));
mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
title(['Real part of the rotated complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
legend('Postive value', 'Negative value');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');

%%%%%%%%%%%%%%%%%%%% Plot imaginary part %%%%%%%%%%%%%%%%%%%%%%
Ynm_rotated_imag = imag(Ynm_rotated);    
figure;
% plot positive values
Ynm_imag_pos = Ynm_rotated_imag.*(Ynm_rotated_imag>=0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_imag_pos);
mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
hold on
% plot negative values
Ynm_imag_neg = Ynm_rotated_imag.*(Ynm_rotated_imag<0);
[Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_imag_neg));
mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
title(['Imaginary part of the rotated complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
legend('Postive value', 'Negative value');
grid on;
axis equal;
xlabel('x');
ylabel('y');
zlabel('z');

end