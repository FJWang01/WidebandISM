
function bn_scatterer_surface = sph_rigid_scatterer_surface(n, a_scaled)
% This function calculates the radial function of the total field on the
% surface of the rigid spherical scatterer of radius a_scaled/k
%
% Formula
% bn(ka) = -1i./(ka).^2./diff_sph_Hankel_2(n, ka) in which ka = a_scaled
%
% Input
% n - order of bn(a_scaled)
% a_scaled - a_scaled/k is the radius of the rigid spherical scatterer
% 
% Note
% n and a_scaled must be of the same size
%
% Output
% bn - bn(a_scaled)
%      size(bn) = size(n) = size(a_scaled)

%% Check if n and a_scaled are of the same size
if ~isequal(size(n), size(a_scaled))
    error('@@ sph_rigid_scatterer_surface: n and a_scaled must be of the same size');
else
    % do nothing
end

%% Check if all values of n are nonnegative integer
validateattributes(n, {'double'}, {'integer', 'nonnegative'});

%% Main
bn_scatterer_surface = -1i./a_scaled.^2./diff_sph_Hankel_2(n, a_scaled);
end