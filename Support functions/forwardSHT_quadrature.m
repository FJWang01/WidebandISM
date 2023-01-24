
function fnm = forwardSHT_quadrature(f, elev, azim, N, weights)
% fnm = forwardSHT_quadrature(f, elev, azim, N, weights)
%
% This function calculates the spherical harmonic coefficients of
% function f sampled at (elev, azim) using quadrature rules. N is the
% truncation order and is a scalar. 
% 
% Inputs: 
% f - function f sampled at Q number of points specified in (elev, azim)
%     size(f) = [Q, T]
% elev - elevation angles of sampling points, in radians
%        size(theta) = [Q, 1]
% azim - azimuth angles of sampling points, in radians
%        size(phi) = [Q, 1]
% N - spherical harmonic truncation order, must be a scalar
% weights - vector of Q weights for the sampling points
%
% Output:
% fnm - spherical harmonic coefficients of function f
%       size(fnm) = [(N+1)^2, T]
%
% Adapted from https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/directSHT.m

%% Check if the number of rows in f is equal to the numbers of elements in elev and azim
if ~isequal(size(f, 1), numel(elev), numel(azim))
    error(['@@ forwardSHT_quadrature: the number of rows of f must be equal ' ...
        'to the numbers of elements in elev and azim']);
else
    % do nothing
end

%% Validate attributes
% elev and azim must be column vectors
validateattributes(elev, {'double'}, {'column'});
validateattributes(azim, {'double'}, {'column'});

validateattributes(N, {'double'}, {'scalar', 'integer', 'nonnegative'});
%% Main
n = repelem((0:N), 2.*(0:N)+1).'; % degree n, (N+1)^2 by 1 column vector
m = (1:(N+1)^2).' - n.^2 - n - 1; % order m, (N+1)^2 by 1 column vector

% size(Ynm_3D) = [numel(n), 1, numel(elev)]
Ynm_3D = sphHarm_mat(n, m, elev, azim);

% size(Ynm_2D) = [numel(n), numel(elev)]
Ynm_2D = squeeze(Ynm_3D); 

Q = size(elev, 1); % number of sampling points

% non-weighted case for uniform/near-uniform arrangements
if nargin<5 || ~exist('weights','var') || isempty(weights)        
    % perform transform
    fnm = (4*pi/Q) * conj(Ynm_2D) * f;
else        
    % perform weighted transform
    fnm = conj(Ynm_2D) * diag(weights) * f;
end

end
