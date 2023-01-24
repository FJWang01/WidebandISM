
function bn_scatt_surface = sph_rigid_scatterer_surface_adapted(n, k, a)
% function bn_scatterer_surface = sph_rigid_scatterer_surface_adapted(n, a_scaled)
%
% This function calculates the radial function for the pressure field on the
% surface of rigid spherical scatterer of radius a
%
% Inputs:
% n - order of bn(ka), must be a column vector
% k - wavenumber, must be a row vector
% a - radius of rigid spherical scatterer, must be a scalar or a column
%     vector
% 
%
% Output:
% bn_scatterer_surface - bn(ka)
%                 size(bn_scatterer_surface) = [numel(n), numel(k), numel(a)]

%% Check the dimensions of inputs
validateattributes(n, {'double'}, {'integer', 'column', 'nonnegative'});
validateattributes(k, {'double'}, {'row'});
validateattributes(a, {'double'}, {'column'});

%% Find unique n values
[unique_n, ~, ic] = unique(n, 'stable');

%% Calculate bn(ka) for unique n
[k_mat, unique_n_mat, a_mat] = meshgrid(k, unique_n, a);

b_scatt_surface_unique_n = sph_rigid_scatterer_surface(unique_n_mat, k_mat .* a_mat);

bn_scatt_surface = b_scatt_surface_unique_n(ic, :, :);
end