
function SH_coeff_new = higher_order_source_in_3D_space_delayed(r, r_s, elev_s, azim_s, v, u, N, time_vec, c, tau)
% SH_coeff_new = higher_order_source_in_3D_space_delayed(r, r_s, elev_s, azim_s, v, u, N, time_vec, c, tau)
%
% This function calculates the spherical harmonic coefficients of the
% signal on the observation sphere due to a source located at (r_s, elev_s,
% azim_s) with directivity patter Yvu.
%
% Inputs:
% r - radius of the observation sphere, in metres, scalar
% r_s - radius of the higher-order source, in metres, scalar
% elev_s - elevation of the higher-order source, in radians, scalar
% azim_s - azimuth of the higher-order source, in radians, scalar
% v - degree of the higher-order source directivity pattern, scalar
% u - order of the higher-order source directivity pattern, scalar
% N - spherical harmonic truncation order, integer scalar
% time_vec - time sample vector, in seconds, must be a row vector
% c - speed of sound, in metres per second, scalar
% tau - delay, in seconds, scalar. The source emits the spherical wave
%       front at t = tau seconds
%
% Output:
% SH_coeff_new - spherical harmonic coefficients of the signal on the
%                observation sphere due to a source of degree v and order u
%                located at spherical coordinate (r_s, elev_s, azim_s)
%                size(SH_coeff)= [(N+1)^2, numel(time_vec)]
%
% Dependencies:
% find_ring_elev_wrt_origin.m
% find_ring_elev_wrt_source.m
% higher_order_source_on_z_axis.m - spherical harmonic coefficients of the
%                                   signal on the observation sphere if the 
%                                   source is located on z-axis
% WignerDM.m - Wigner D-matrix to rotate the signal on the observation sphere

%% Validate attributes
validateattributes(r, {'double'}, {'scalar'});
validateattributes(r_s, {'double'}, {'scalar'});
validateattributes(elev_s, {'double'}, {'scalar'});
validateattributes(azim_s, {'double'}, {'scalar'});
validateattributes(v, {'double'}, {'scalar', 'integer', 'nonnegative'}); % degree v
validateattributes(u, {'double'}, {'scalar', 'integer'}); % order u
validateattributes(N, {'double'}, {'scalar', 'integer', 'nonnegative'}); % truncation order N
validateattributes(c, {'double'}, {'scalar'});
validateattributes(time_vec, {'double'}, {'row'}); % must be row vector
validateattributes(tau, {'double'}, {'scalar', 'nonnegative'});

%% Find SH coefficients of the signal on the observation sphere when the source is on z-axis
SH_coeff = higher_order_source_on_z_axis_delayed(r, r_s, v, u, N, time_vec, c, tau);

%% Rotate the observed signal first by elev_s about y-axis, then by azim_s about z-axis
Wigner_D_mat = WignerDM(N, azim_s, elev_s, 0);
SH_coeff_new = Wigner_D_mat*SH_coeff;
end