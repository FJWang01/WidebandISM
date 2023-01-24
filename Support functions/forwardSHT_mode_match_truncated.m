
function [n, m, nm_mask, fnm] = forwardSHT_mode_match_truncated(f, elev, azim, N)
% [n, m, nm_mask, fnm] = forwardSHT_mode_match_truncated(f, elev, azim, N)
%
% This function calculates the spherical harmonic coefficients of signal f
% sampled at angular directions (elev, azim). The spherical harmonic
% truncation orders are in row vector N. numel(N) = size(f, 2).
%
% Inputs:
% f - signal f sampled at angular directions (elev, azim)
% elev - elevation angles of sampling points, must be a column vector
% azim - azimuth angles of sampling points, must be a column vector
% N - spherical harmonic truncation order, must be a row vector
%     numel(N) = size(f, 2)
%
% Notes:
% 1. elev and azim must be of the same size and must be column vectors
% 2. size(f, 1) = numel(elev) = numel(azim)
% 3. N must be a row vector and numel(N) = size(f, 2)
%
% Outputs:
% n - degrees of SH coefficients correspond to fnm
% m - orders of SH coefficients correspond to fnm
% nm_mask - mask to mask out unused n and m
%           size(n) = size(m) = size(nm_mask) = [(max(N)+1)^2, size(f, 2)]
% fnm - spherical harmonic coefficients of signal f
%       size(fnm) = [(max(N)+1)^2, size(f, 2)]

%% Check dimenstions of inputs
% Check if the number of rows of f corresponds to the numbers of elements
% in elev and azim
if ~isequal(size(f, 1), numel(elev), numel(azim))
    error(['@@ forwardSHT_mode_match_truncated: the number of rows of f ' ...
        'must be equal to the numbers of elements in elev and azim']);
else
    % do nothing
end

% Check if the number of columns of f corresponds to the number of elements
% in N
if ~isequal(size(f, 2), numel(N))
    error(['@@ forwardSHT_mode_match_truncated: the number of columns of f ' ...
        'must be equal to the number of elements in N']);
else
    % do nothing
end

%% Validate attributes
validateattributes(elev, {'double'}, {'column'});
validateattributes(azim, {'double'}, {'column'});

% SH truncation orders N must be a row vector
validateattributes(N, {'double'}, {'integer', 'row', 'nonnegative'}); 

%% Generate degree-order pair correspond to max(N)
N_max = max(N);

% size(n_max) = size(m_max) = [(max(N)+1)^2, 1]
[n_max, m_max, ~] = getDegreeOrderPairs(N_max);

%% Get n, m and nm_mask
[n, m, nm_mask] = getDegreeOrderPairs(N);

%% Generate mode matching matrix
% size(Ynm) = [(max(N)+1)^2, 1, numel(elev)]
Ynm = sphHarm_mat(n_max, m_max, elev, azim); 

% size(Ynm_2D) = [(max(N)+1)^2, numel(elev)]
Ynm_2D = squeeze(Ynm);

% size(Ynm_2D_final) = [numel(elev), (max(N)+1)^2]
Ynm_2D_final = Ynm_2D.';

fnm = zeros((N_max + 1)^2, size(f, 2));
for k_idx = 1:size(f, 2)
    truncation_lim = (N(k_idx)+1)^2;
    Ynm_2D_final_truncated = Ynm_2D_final(:, 1:truncation_lim);
    fnm(1:truncation_lim, k_idx) = pinv(Ynm_2D_final_truncated) * f(:, k_idx);
end
end