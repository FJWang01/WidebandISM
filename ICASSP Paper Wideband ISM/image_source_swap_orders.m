
function fnm_swapped = image_source_swap_orders(fnm, swap_mat, scalings)
% This function swaps the spherical harmonic coefficients of the opposite
% order
%
% Inputs
% fnm - SH coefficients of directivities of spherical wave fronts emitted
%       by the original source
%       size(fnm, 2) = numel(tau) where tau are the time instances when
%       spherical wave fronts are emitted
% swap_mat - if swapping spherical harmonic coefficients of the opposite
%            order is needed. Logical 1 means yes, Logical zero means no
%            numel(swap_mat) = number of image sources
% scalings - -1 or 1 scaling due to reflections
%            size(scalings) = [size(fnm, 1), numel(swap_mat)]
%
% Output
% fnm_swapped - SH coefficients of directivities of spherical wave fronts
%               emitted by the image sources 
%               image sources are the second dimension, tau is the third
%               dimension
%               size(fnm_swapped) = [size(fnm, 1), numel(swap_mat), size(fnm, 2)]

%% Check inputs
validateattributes(swap_mat, {'logical'},{'row'});

if ~isequal(size(fnm, 1), size(scalings, 1)) || ~isequal(numel(swap_mat), size(scalings, 2))
    error('@@ image_source_swap_orders: the dimension of +1/-1 scaling is incorrect');
else
    % do nothing
end

%% Main
% tau becomes the 3rd dimension
fnm_reshaped = reshape(fnm, size(fnm, 1), 1, size(fnm, 2));

% Apply scaling (1 or -1) due to reflection 
% size(fnm_3D) = [size(fnm, 1), numel(swap_mat), size(fnm, 2)]
% size(scalings) = [size(fnm, 1), numel(swap_mat)]
fnm_3D = repmat(fnm_reshaped, 1, numel(swap_mat)); 
fnm_scaled = fnm_3D .* scalings; 

% Swap SH coefficients of the opposite order according to swap_mat
% Initialise fnm_swapped
fnm_swapped = fnm_scaled; 
N = sqrt(size(fnm, 1))-1;
for n = 0:N
    fnm_swapped(n^2+1:n^2+2*n+1, swap_mat, :) = flip(fnm_scaled(n^2+1:n^2+2*n+1, swap_mat, :), 1);
end
end