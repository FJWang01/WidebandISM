
function bn = sph_rigid_scatterer_adapted(n, k, r, a)
% bn = sph_rigid_scatterer_adapted(n, k, r, a)
%
% This function calculates the radial function bn(kr)
% 
% Inputs:
% n - order of bn, must be a column vector
% k - wavenumber, must be a row vector
% r - radius, must be a column vector
% a - radius of the spherical rigid scatterer, must be a scalar
%
% Output:
% bn - size(bn) = [numel(n), numel(k), numel(r)]

%% Validate attributes of inputs
validateattributes(n, {'double'}, {'integer', 'column', 'nonnegative'});
validateattributes(k, {'double'}, {'row'});
validateattributes(r, {'double'}, {'column'});
validateattributes(a, {'double'}, {'scalar'});

%% Find unique n values
[unique_n, ~, ic] = unique(n, 'stable');

%% Calculate bn(kr) for unique n
[k_mat, unique_n_mat, r_mat] = meshgrid(k, unique_n, r);

b_unique_n = sph_rigid_scatterer(unique_n_mat, k_mat .* r_mat, k_mat * a);

bn = b_unique_n(ic, :, :);

end
