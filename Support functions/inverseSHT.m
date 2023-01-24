
function f = inverseSHT(fnm, elev, azim, n, m)
% f = inverseSHT(fnm, elev, azim, n, m)
% 
% This function calculates the values of function f at sample locations
% (elev, azim) based on the spherical harmonic coefficients
%
% Inputs:
% fnm - spherical harmonic coefficients of function f
% elev - elevation angles of sample locations, must be a column vector
% azim - azimuth angles of sample locations, must be a column vector
% n - degree of spherical harmonic functions, must be a column vector
% m - order of spherical harmonic functions, must be a column vector
%
% Notes:
% 1. fnm, n and m must have the same number of rows
% 2. n and m nust be of the same size and must be column vectors
% 3. elev and azim must be of the same size and must be column vectors
%
% Output:
% f - the values of function f at sample locations (elev, azim)
%     size(f) = [numel(elev), size(fnm, 2)]

%% Check if the number of rows of fnm is equal to the numbers of elements in n and m
if ~isequal(size(fnm, 1), numel(n), numel(m))
    error(['@@ inverseSHT: the number of rows of fnm must be equal to the ' ...
        'numbers of elemements in n and m']);
else
    % do nothing
end

%% Validate attributes
validateattributes(n, {'double'}, {'integer', 'nonnegative', 'column'});
validateattributes(m, {'double'}, {'integer', 'column'});

validateattributes(elev, {'double'}, {'column'});
validateattributes(azim, {'double'}, {'column'});

%% Calculate Ynm
% size(Ynm_3D) = [numel(n), 1, numel(elev)]
Ynm_3D = sphHarm_mat(n, m, elev, azim);

% size(Ynm_2D) = [numel(n), numel(elev)]
Ynm_2D = squeeze(Ynm_3D);

% size(Ynm_2D_final) = [numel(elev), numel(n)]
Ynm_2D_final = Ynm_2D.';

%% Calculate f
% size(f) = [numel(elev), size(SH_coeff, 2)]
f = Ynm_2D_final * fnm;
end