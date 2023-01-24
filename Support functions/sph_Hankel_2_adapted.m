
function h_Han_2 = sph_Hankel_2_adapted(n, k, r)
% h_Han_2 = sph_Hankel_2_adapted(n, k, r)
%
% This function calculates the spherical Hankel function of the 2nd kind
% hn^(2)(kr)
% 
% Inputs:
% n - order of the hn^(2)(kr), must be a column vector
% k - wavenumber, must be a row vector
% r - radius, must be a column vector
%
% Output:
% h_Han_2 - values of hn^(2)(kr)
%           size(h_Han_2) = [numel(n), numel(k), numel(r)]

%% Validate attributes of inputs
validateattributes(n, {'double'}, {'integer', 'column', 'nonnegative'});
validateattributes(k, {'double'}, {'row'});
validateattributes(r, {'double'}, {'column'});

%% Find unique n values
[unique_n, ~, ic] = unique(n, 'stable');

%% Calculate jn(kr) for unique n
[k_mat, unique_n_mat, r_mat] = meshgrid(k, unique_n, r);

h_Han_2_unique_n = sph_Hankel_2(unique_n_mat, k_mat .* r_mat);

h_Han_2 = h_Han_2_unique_n(ic, :, :);

end
