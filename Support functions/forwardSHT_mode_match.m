
function fnm = forwardSHT_mode_match(f, elev, azim, N)
% fnm = forwardSHT_mode_match(f, elev, azim, N)
%
% This function calculates the spherical harmonic coefficients of signal f
% sampled at angular directions (elev, azim). The spherical harmonic
% truncation order is N. N is a scalar. 
%
% Inputs:
% f - signal f sampled at angular directions (elev, azim)
% elev - elevation angles of sampling points, must be a column vector
% azim - azimuth angles of sampling points, must be a column vector
% N - spherical harmonic truncation order, must be a scalar 
%
% Notes:
% 1. elev and azim must be of the same size and must be column vectors
% 2. size(f, 1) = numel(elev) = numel(azim)
% 3. N must be a scalar
%
% Output:
% fnm - spherical harmonic coefficients of signal f
%       size(fnm) = [(N+1)^2, size(f, 2)]

%% Check if the number of rows of f is equal to the numbers of elements in elev and azim
if ~isequal(size(f, 1), numel(elev), numel(azim))
    error(['@@ forwardSHT_mode_match: the number of row of f must be equal ' ...
        'to the numbers of elements in elev and azim']);
else
    % do nothing
end

%% Validate attributes
validateattributes(elev, {'double'}, {'column'});
validateattributes(azim, {'double'}, {'column'});

validateattributes(N, {'double'}, {'scalar', 'integer', 'nonnegative'});

%% Generate degree-order pair
n = (repelem((0:N), 2*(0:N)+1)).'; % degrees, (N+1)^2 by 1 column vector
m = (1:(N+1)^2).' - n.^2 - n - 1; % orders, (N+1)^2 by 1 column vector

%% Generate mode matching matrix
Ynm = sphHarm_mat(n, m, elev, azim); % size(Ynm) = [(N+1)^2, 1, numel(elev)]
Ynm_2D = squeeze(Ynm); % size(Ynm_2D) = [(N+1)^2, numel(elev)]
Ynm_2D_final = Ynm_2D.'; % size(Ynm_2D_final) = [numel(elev), (N+1)^2]

%% Solve for fnm
fnm = pinv(Ynm_2D_final) * f;
end