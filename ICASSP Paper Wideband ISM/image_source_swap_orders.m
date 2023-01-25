
function fnm_final = image_source_swap_orders(fnm, swap_mat, scalings, refl_coeff)
% fnm_final = image_source_swap_orders(fnm, swap_mat, scalings, refl_coeff)
%
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
% refl_coeff - wall reflection coefficients of the image sources
%              size(refl_coeff) = [1, numel(swap_mat)]
%
% Output
% fnm_final - SH coefficients of directivities of spherical wave fronts
%             emitted by the image sources 
%             image sources are the second dimension, tau is the third
%             dimension
%             size(fnm_final) = [size(fnm, 1), numel(swap_mat), size(fnm, 2)]

%% Check inputs
validateattributes(swap_mat, {'logical'}, {'row'});
validateattributes(refl_coeff, {'double'}, {'row'});

if ~isequal(size(fnm, 1), size(scalings, 1)) || ~isequal(numel(swap_mat), size(scalings, 2))
    error('@@ image_source_swap_orders: the dimension of +1/-1 scaling is incorrect');
else
    % do nothing
end

if ~isequal(numel(swap_mat), numel(refl_coeff))
    error('@@ image_source_swap_orders: the dimension of swap_mat and refl_coeff must be equal');
else
    % do nothing
end

%% Main
% size(fnm_reshaped) = [size(fnm, 1), 1, size(fnm, 2)]
% tau is the third dimension
fnm_reshaped = reshape(fnm, size(fnm, 1), 1, size(fnm, 2));

% size(fnm_3D) = [size(fnm, 1), numel(swap_mat), size(fnm, 2)]
fnm_3D = repmat(fnm_reshaped, 1, numel(swap_mat));

% Apply scaling (1 or -1) due to reflection 
% size(scalings) = [size(fnm, 1), numel(swap_mat)]
% size(fnm_scaled) = [size(fnm, 1), numel(swap_mat), size(fnm, 2)]
fnm_scaled = fnm_3D .* scalings; 

% Swap SH coefficients of the opposite order according to swap_mat
% Initialise fnm_swapped
fnm_swapped = fnm_scaled;

% Get the max SH truncation order of fnm
N = sqrt(size(fnm, 1))-1;

% Go through each truncation order and swap opposite order as desired
for n = 0:N
    fnm_swapped(n^2+1:n^2+2*n+1, swap_mat, :) = flip(fnm_scaled(n^2+1:n^2+2*n+1, swap_mat, :), 1);
end

% Apply wall reflection coefficients
fnm_final = refl_coeff .* fnm_swapped; 
end