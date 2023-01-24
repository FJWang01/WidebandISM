
function cos_theta_dash = find_ring_elev_wrt_source(r, r_s, c, time_vec)
% [cos_theta_dash] = find_ring_elev_wrt_source(r, r_s, c, time_vec)
%
% This function finds the elevation angle of the ring which is the
% intersection of the observation sphere and the expanding spherical 
% surface with respect to the origin of the source
%
% Inputs:
% r - radius of the observation sphere, in metres, must be a scalar
% r_s - source radius, in metres, must be a scalar
% c - speed of sound, in metres per second, scalar
% time_vec - time sample vector, in seconds, must be a row vector
% 
% Output:
% cos_theta_dash - cosine of the elevation angle of the ring w.r.t. the 
%                  origin of the source
%                  size(cos_theta_dash) = [1, numel(time_vec)]

%% Validate attributes
% Check if r, r_s and c are scalars
validateattributes(r, {'double'}, {'scalar', 'nonnegative'});
validateattributes(r_s, {'double'}, {'scalar', 'nonnegative'});
validateattributes(c, {'double'}, {'scalar'});

% Check if time_vec is a row vector
validateattributes(time_vec, {'double'}, {'row'});

%% Mains
% range of time during which intersection occurs
range = (time_vec >= abs(r-r_s)/c) & (time_vec <= (r+r_s)/c); 

% alpha = pi - theta_dash
cos_alpha = zeros(1, numel(time_vec)); % shape dimension
cos_alpha(range) = (c^2*time_vec(range).^2 + r_s^2 - r^2)./(2*c*time_vec(range) * r_s);
cos_theta_dash = -cos_alpha; % cos(pi - theta_dash) = -cos(theta_dash)

end