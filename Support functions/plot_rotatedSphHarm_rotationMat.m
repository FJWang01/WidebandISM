
function [] = plot_rotatedSphHarm_rotationMat(n, m, alpha, beta, gamma)
% This function plots the rotated spherical harmonic function of degree n
% and order m using the rotation matrix

% Input
% n - degree of the spherical harmonic function, scalar
% m - order of the spherical harmonic function, scalar
% alpha, beta, gamma - rotation angles in radian

% Output
% None

% Rotations follow zyz convention
% First, rotate gamma radian about z-axis
% Next, rotate beta radian about y-axis
% Then, rotate alpha radian about z-axis
%% Plotting grid
[X, Y, Z] = sphere(128);
[azim, elev, ~] = cart2sph(X, Y, Z);
elev = pi/2 - elev; 
%% Original spherical harmonic function
Ynm = sphHarm(n, m, elev, azim);
%% Get rotation matrix
R = RZRYRZrad(alpha, beta, gamma);
%% Rotate real part
Ynm_real = real(Ynm);
% real positive
Ynm_real_pos = Ynm_real.*(Ynm_real>=0);
[x_real_pos, y_real_pos, z_real_pos] = sph2cart(azim, pi/2 - elev, Ynm_real_pos);

xyz_rotated_real_pos = R*[x_real_pos(:) y_real_pos(:) z_real_pos(:)]';

x_real_pos(:) = xyz_rotated_real_pos(1,:)'; % refill x_real_pos matrix
y_real_pos(:) = xyz_rotated_real_pos(2,:)'; % refill y_real_pos matrix
z_real_pos(:) = xyz_rotated_real_pos(3,:)'; % refill z_real_pos matrix

% real negative
Ynm_real_neg = Ynm_real.*(Ynm_real<0);
[x_real_neg, y_real_neg, z_real_neg] = sph2cart(azim, pi/2 - elev, abs(Ynm_real_neg));

xyz_rotated_real_neg = R*[x_real_neg(:) y_real_neg(:) z_real_neg(:)]';

x_real_neg(:) = xyz_rotated_real_neg(1,:)'; % refill x_real_neg matrix
y_real_neg(:) = xyz_rotated_real_neg(2,:)'; % refill y_real_neg matrix
z_real_neg(:) = xyz_rotated_real_neg(3,:)'; % refill z_real_neg matrix

figure;
mesh(x_real_pos, y_real_pos, z_real_pos, 'EdgeColor', 'c');
hold on
mesh(x_real_neg, y_real_neg, z_real_neg, 'EdgeColor', 'm');
legend('Postive value', 'Negative value');
title(['Real part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;

%% Rotate imaginary part
Ynm_imag = imag(Ynm);
% imaginary positive
Ynm_imag_pos = Ynm_imag.*(Ynm_imag>=0);
[x_imag_pos, y_imag_pos, z_imag_pos] = sph2cart(azim, pi/2 - elev, Ynm_imag_pos);

xyz_rotated_imag_pos = R*[x_imag_pos(:) y_imag_pos(:) z_imag_pos(:)]';

x_imag_pos(:) = xyz_rotated_imag_pos(1,:)'; % refill x_imag_pos matrix
y_imag_pos(:) = xyz_rotated_imag_pos(2,:)'; % refill y_imag_pos matrix
z_imag_pos(:) = xyz_rotated_imag_pos(3,:)'; % refill z_imag_pos matrix

% imaginary negative
Ynm_imag_neg = Ynm_imag.*(Ynm_imag<0);
[x_imag_neg, y_imag_neg, z_imag_neg] = sph2cart(azim, pi/2 - elev, abs(Ynm_imag_neg));

xyz_rotated_imag_neg = R*[x_imag_neg(:) y_imag_neg(:) z_imag_neg(:)]';

x_imag_neg(:) = xyz_rotated_imag_neg(1,:)'; % refill x_imag_neg matrix
y_imag_neg(:) = xyz_rotated_imag_neg(2,:)'; % refill y_imag_neg matrix
z_imag_neg(:) = xyz_rotated_imag_neg(3,:)'; % refill z_imag_neg matrix

figure;
mesh(x_imag_pos, y_imag_pos, z_imag_pos, 'EdgeColor', 'c');
hold on
mesh(x_imag_neg, y_imag_neg, z_imag_neg, 'EdgeColor', 'm');
legend('Postive value', 'Negative value');
title(['Imaginary part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
end

