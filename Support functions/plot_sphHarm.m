
function Ynm = plot_sphHarm(n, m, type)
% This function plots the spherical harmonic function of degree n and order m
% Input
% n - degree, scalar
% m - order, scalar
% type - 'R' for plotting real spherical harmonics or 'C' for plotting complex
% spherical harmonics

%% Generate sampling points on unit sphere
[X, Y, Z] = sphere(128);
[azim, elev, ~] = cart2sph(X, Y, Z);
elev = pi/2 - elev; 
%% Convert to Cartesian coordinate and plot as desired
if type == 'C' % complex spherical harmonics
    Ynm = sphHarm(n, m, elev, azim);
    
    %%%%%%%%%%%%%%%%%%%% Plot real part %%%%%%%%%%%%%%%%%%%%%%
    Ynm_real = real(Ynm);    
    figure;
    % plot positive values
    Ynm_real_pos = Ynm_real.*(Ynm_real>=0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_real_pos);
    mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
    hold on
    % plot negative values
    Ynm_real_neg = Ynm_real.*(Ynm_real<0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_real_neg));
    mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
    title(['Real part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
    legend('Postive value', 'Negative value');
    grid on;
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
    %%%%%%%%%%%%%%%%%%%% Plot imaginary part %%%%%%%%%%%%%%%%%%%%%%
    Ynm_imag = imag(Ynm);    
    figure;
    % plot positive values
    Ynm_imag_pos = Ynm_imag.*(Ynm_imag>=0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_imag_pos);
    mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
    hold on
    % plot negative values
    Ynm_imag_neg = Ynm_imag.*(Ynm_imag<0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_imag_neg));
    mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
    title(['Imaginary part of complex spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
    legend('Postive value', 'Negative value');
    grid on;
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    
elseif type == 'R' % real spherical harmonics
    Ynm = realSphHarm(n, m, elev, azim); 
    figure;
    % plot positive values
    Ynm_pos = Ynm.*(Ynm>=0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, Ynm_pos);
    mesh(Xm, Ym, Zm, 'EdgeColor', 'c');
    hold on
    % plot negative values
    Ynm_neg = Ynm.*(Ynm<0);
    [Xm, Ym, Zm] = sph2cart(azim, pi/2 - elev, abs(Ynm_neg));
    mesh(Xm, Ym, Zm, 'EdgeColor', 'm');
    title(['Real spherical harmonic function of degree ', num2str(n), ' and order ', num2str(m)]);
    legend('Postive value', 'Negative value');
    grid on;
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
else
    error('Please select real (R) or complex (C) spherical harmonics');
end
end