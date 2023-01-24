
function bn = sph_rigid_scatterer(n, r_scaled, a_scaled)
% This function calculates the radial function of a rigid spherical
% scatterer bn(r_scaled)
%
% Input
% n - order of bn(r_scaled)
% r_scaled - argument, r_scaled/k is the radius of interest
% a_scaled - a_scaled/k is the radius of the rigid spherical scatterer
% 
% Notes
% 1. n, r_scaled and a_scaled must be of the same size
% 2. r_scaled must be greater than or equal to a_scaled
%
% Output
% bn - bn(r_scaled)
%      size(bn) = size(n) = size(r_scaled) = size(a_scaled)

%% Check if n, r_scaled and a_scaled are of the same size
if ~isequal(size(n), size(r_scaled), size(a_scaled))
    error('@@ sph_rigid_scatterer: n, r_scaled and a_scaled must be of the same size');
else
    % do nothing
end

%% Check if all values in r_scaled are greater than or equal to corresponding values in a_scaled
if any(r_scaled < a_scaled, 'all')
    error('@@ sph_rigid_scatterer: r_scaled must be greater than or equal to a_scaled');
else
    % do nothing
end

%% Validate attributes
validateattributes(n, {'double'}, {'integer', 'nonnegative'});

%% Calculate bn(r_scaled)
bn = sph_Bessel_1(n, r_scaled) - ...
    diff_sph_Bessel_1(n, a_scaled)./diff_sph_Hankel_2(n, a_scaled).*sph_Hankel_2(n, r_scaled);
end