
function y_bes = sph_Bessel_2_adapted(n, k, r)
% y_bes = sph_Bessel_2_adapted(n, k, r)
%
% This function calculates the spherical Bessel function of the 2nd kind
% yn(kr)
% 
% Inputs:
% n - order of the yn(kr), must be a column vector
% k - wavenumber, must be a row vector
% r - radius, must be a column vector
%
% Output:
% y_bes - values of yn(kr)
%         size(y_bes) = [numel(n), numel(k), numel(r)]

%% Validate attributes of inputs
validateattributes(n, {'double'}, {'integer', 'column', 'nonnegative'});
validateattributes(k, {'double'}, {'row'});
validateattributes(r, {'double'}, {'column'});

%% Find unique n values
[unique_n, ~, ic] = unique(n, 'stable');

%% Calculate jn(kr) for unique n
[k_mat, unique_n_mat, r_mat] = meshgrid(k, unique_n, r);

y_bes_unique_n = sph_Bessel_2(unique_n_mat, k_mat .* r_mat);

y_bes = y_bes_unique_n(ic, :, :);

end
