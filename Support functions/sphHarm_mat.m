
function Ynm = sphHarm_mat(n, m, theta, phi)
% This function calculates the spherical harmonic functions of degree n and
% order m at locations specified by theta and phi
%
% Input
% n - degree
% m - order
% theta - elevation, in radian
% phi - azimuth, in radian
% 
% Note:
% 1. n and m must be of the same size
% 2. theta and phi must be of the same size and must be column vectors
%
% Output
% Ynm - size(Ynm) = [size(n, 1), size(n, 2), numel(theta)]

%% Check if n and m are of the same size
if ~isequal(size(n), size(m))
    error('@@ sphHarm_mat: n and m must be of the same size');
else
    % do nothing
end

%% Check if theta and phi are of the same size
if ~isequal(size(theta), size(phi))
    error('@@ sphHarm_mat: theta and phi must be of the same size');
else
    % do nothing
end

%% Check if the values of n and m are valid
validateattributes(n, {'double'}, {'integer', 'nonnegative'});
validateattributes(m, {'double'}, {'integer'});

%% Check if theta and phi are column vectors
validateattributes(theta, {'double'}, {'column'});
validateattributes(phi, {'double'}, {'column'});

%% Calculate Pnm(cos(theta))
Pnm = assoc_Legendre_poly(n, m, cos(theta.')); % size(Pnm) = [size(n, 1), size(n, 2), numel(theta)]

%% Calculate exp(i*m*phi)
m_3D = repmat(m, 1, 1, numel(phi)); % size(m_3D) = [size(m, 1), size(m, 2), numel(phi)]
phi_3D = reshape(phi, 1, 1, numel(phi)); % size(phi_3D) = [1, 1, numel(phi)]
exp_term = exp(1i.*m_3D.*phi_3D); % size(exp_term) = [size(m, 1), size(m, 2), numel(phi)]

%% Calculate the normalisation factor
Anm = sqrt((2*n+1).*factorial(n-m)./factorial(n+m)/4/pi); % size(Anm) = size(n)

%% Calculate Ynm
Ynm = Anm.*Pnm.*exp_term; % size(Ynm) = [size(n, 1), size(n, 2), numel(theta)]
end