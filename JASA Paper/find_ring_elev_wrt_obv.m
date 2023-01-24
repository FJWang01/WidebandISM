function cos_theta = find_ring_elev_wrt_obv(r, r_s, c, time_vec)
% [cos_theta] = find_ring_elev_wrt_obv(r, r_s, c, time_vec)
%
% This function finds the elevation angle of the ring which is the
% intersection of the observation sphere and the expanding spherical 
% surface with respect to the origin of the observation sphere
%
% Inputs:
% r - radius of the observation sphere, in metres, must be a scalar
% r_s - source radius, in metres, must be a scalar
% c - speed of sound, in metres per second, scalar
% time_vec - time sample vector, in seconds, must be a row vector
%
%
% Output:
% cos_theta - cosine of the elevation angle of the ring w.r.t. the 
%             origin of the observation sphere
%             size(cos_theta) = [1, numel(time_vec)]

%% Validate attributes
% Check if r, r_s and c are scalars
validateattributes(r, {'double'}, {'scalar', 'nonnegative'});
validateattributes(r_s, {'double'}, {'scalar', 'nonnegative'});
validateattributes(c, {'double'}, {'scalar'});

% Check if time_vec is a row vector
validateattributes(time_vec, {'double'}, {'row'});

%% Mains
% Range of time during which intersection occurs
range = (time_vec >= abs(r-r_s)/c) & (time_vec <= (r+r_s)/c); 

% Find cos_theta, using the law of cosines
cos_theta = zeros(1, numel(time_vec)); % shape dimension 
cos_theta(range) = (r^2 + r_s^2 - c^2*time_vec(range).^2)/(2*r*r_s);

end