
function Ynm_rotated = getRotatedSphHarm(n, m, alpha, beta, gamma, theta, phi)
% This function calculates the rotated sphericla harmonic function of
% degree n and order m at locations specified by theta and phi

% Input
% n - degree of the spherical harmonic function
% m - order of the spherical harmonic function
% (theta, phi) - locations of points at which the rotated spherical
%                harmonic function are calculated
% alpha, beta, gamma - rotation angles, in radian

% Output
% Ynm_rotated - the values of the rotated spherical harmonic function of
%               degree n and order m at locations specified by (theta, phi)
% size(Ynm_rotated) = size(theta) = size(phi)

% Note:
% 1. Rotation follows zyz convention
%     First, rotate gamma by z-axis
%     Next, rotate beta by y-axis
%     Then, rotate alpha by z-axis
% 2. theta and phi must be of the same size

%% Check if theta and phi are of the same size
if ~isequal(size(theta),size(phi))
	error('@@ getRotatedSphHarm: theta and phi must be of the same size');
else 
    % do nothing
end
%% Obtain new spherical harmonic coefficients
% Wigner D-matrix
N = n; % spherical harmonic truncation order for generating Wigner D-matrix
DM = WignerDM(N, alpha, beta, gamma);

% Old spherical harmonic coefficients to be multiplied with the Wigner D-matrix
old_fnm = zeros((N+1)^2, 1);
old_fnm(n^2+n+1+m) = 1;

% Get new spherical harmonic coefficients
new_fnm = DM * old_fnm;
new_fnm_3D = reshape(new_fnm, 1, 1, numel(new_fnm));

%% Calculate the values of the rotated spherical harmonic function at (elev, azim)
new_m = -n:n; % rotation does not mix degrees, it only mixes orders
Ynm_at_n = zeros(size(theta, 1), size(theta, 2), (N+1)^2); 

for idx = 1:numel(new_m) % go through each order at the desired degree n
    Ynm = sphHarm(n, new_m(idx), theta, phi); % size(Ynm) = size(theta) = size(phi)
    Ynm_at_n(:, :, (N)^2 + idx) = Ynm;
end

Ynm_rotated_ind_mode = Ynm_at_n .* new_fnm_3D;

Ynm_rotated = sum(Ynm_rotated_ind_mode, 3);

end