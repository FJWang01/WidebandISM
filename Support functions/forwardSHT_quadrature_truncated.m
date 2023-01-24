
function [n, m, nm_mask, fnm] = forwardSHT_quadrature_truncated(f, elev, azim, N, weights)
% [n, m, nm_mask, fnm] = forwardSHT_quadrature_truncated(f, elev, azim, N, weights)
%
% This function calculates the spherical harmonic coefficients of
% function f sampled at (elev, azim) using quadrature rules. The truncation
% orders are in row vector N. numel(N) = size(f, 2).
% 
% Inputs: 
% f - function f sampled at Q number of points specified in (elev, azim)
%     size(f) = [Q, T]
% elev - elevation angles of sampling points, in radians
%        size(theta) = [Q, 1]
% azim - azimuth angles of sampling points, in radians
%        size(phi) = [Q, 1]
% N - spherical harmonic truncation order, must be a row vector
%     size(N) = [1, T]
% weights - vector of Q weights for the sampling points
%
% Outputs:
% n - degrees of SH coefficients correspond to fnm
% m - orders of SH coefficients correspond to fnm
% nm_mask - mask to mask out unused n and m
% fnm - spherical harmonic coefficients of function f
%       size(fnm) = [(max(N)+1)^2, T]
%
% Adapted from https://github.com/polarch/Spherical-Harmonic-Transform/blob/master/directSHT.m

%% Check the dimensions of inputs
% Check if the number of rows of f is equal to the numbers of elements in 
% elev and azim
if ~isequal(size(f, 1), numel(elev), numel(azim))
    error(['@@ forwardSHT_quadrature_truncated: the number of rows of f ' ...
        'must be equal to the numbers of elements in elev and azim']);
else
    % do nothing
end

% Check if the number of columns of f is equal to the number of elements in
% N
if ~isequal(size(f, 2), numel(N))
    error(['@@ forwardSHT_quadrature_truncated: the number of columns of f ' ...
        'must be equal to the number of elements in N']);
else
    % do nothing
end

%% Validate attributes
% elev and azim must be column vectors
validateattributes(elev, {'double'}, {'column'});
validateattributes(azim, {'double'}, {'column'});

% N must be a row vector
validateattributes(N, {'double'}, {'integer', 'row', 'nonnegative'});

%% Generate degree-order pairs correspond to max(N)
N_max = max(N);
[n_max, m_max, ~] = getDegreeOrderPairs(N_max);

%% Get n, m and nm_mask
[n, m, nm_mask] = getDegreeOrderPairs(N);

%% Calculate fnm
% size(Ynm_3D) = [numel(n_max), 1, numel(elev)]
Ynm_3D = sphHarm_mat(n_max, m_max, elev, azim);

% size(Ynm_2D) = [numel(n), numel(elev)]
Ynm_2D = squeeze(Ynm_3D); 

Q = size(elev, 1); % number of sampling points

fnm = zeros((N_max+1)^2, size(f, 2));
for k_idx = 1:size(f, 2)
    truncation_lim = (N(k_idx)+1)^2;
    % non-weighted case for uniform/near-uniform arrangements
    if nargin<5 || ~exist('weights','var') || isempty(weights)        
        % perform transform
        fnm(1:truncation_lim, k_idx) = (4*pi/Q) * conj(Ynm_2D(1:truncation_lim, :)) * f(:, k_idx);
    else        
        % perform weighted transform
        fnm(1:truncation_lim, k_idx) = conj(Ynm_2D(1:truncation_lim, :)) * diag(weights) * f(:, k_idx);
    end
end

end
