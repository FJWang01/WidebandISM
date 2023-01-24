
function j_bes = sph_Bessel_1_adapted(n, k, r)
% j_bes = sph_Bessel_1_adapted(n, k, r)
%
% This function calculates the spherical Bessel function of the 1st kind
% jn(kr)
% 
% Inputs:
% n - order of the jn(kr), must be a column vector
% k - wavenumber, must be a row vector
% r - radius, must be a column vector
%
% Output:
% j_bes - values of jn(kr)
%         size(j_bes) = [numel(n), numel(k), numel(r)]

%% Validate attributes of inputs
validateattributes(n, {'double'}, {'integer', 'column', 'nonnegative'});
validateattributes(k, {'double'}, {'row'});
validateattributes(r, {'double'}, {'column'});

%% Find unique n values
[unique_n, ~, ic] = unique(n, 'stable');

%% Calculate jn(kr) for unique n
[k_mat, unique_n_mat, r_mat] = meshgrid(k, unique_n, r);

j_bes_unique_n = sph_Bessel_1(unique_n_mat, k_mat .* r_mat);

j_bes = j_bes_unique_n(ic, :, :);

end
