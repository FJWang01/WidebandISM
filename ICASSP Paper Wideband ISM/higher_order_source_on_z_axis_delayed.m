
function SH_coeff = higher_order_source_on_z_axis_delayed(r, r_s, v, u, N, time_vec, c, tau)
% SH_coeff = higher_order_source_on_z_axis_delayed(r, r_s, v, u, N, time_vec, c, tau)
%
% This function calculates the spherical harmonic coefficients of the
% observed signal on the observation sphere due to a source of degree v and 
% order u located on the positive z-axis
%
% Inputs:
% r - radius of the observation sphere, in metres, scalar
% r_s - source radius, in metres, scalar
% v - degree of the source, scalar
% u - order of the source, scalar
% N - spherical harmonic truncation order of the observed signal, scalar
% time_vec - time samples, in seconds, must be a row vector
% c - speed of sound, in metres per second, scalar
% tau - delay, in seconds, scalar. The source emits the spherical wave
%       front at t = tau seconds
%
% Output:
% SH_coeff - spherical harmonic coefficients of the observed signal on the 
%            observation sphere due to a source of degree v and order u 
%            located on the positive z-axis
%            size(SH_coeff) = [(N+1)^2, numel(time_vec)]
%

%% Validate attributes of inputs
validateattributes(r, {'double'}, {'scalar'});
validateattributes(r_s, {'double'}, {'scalar'});
validateattributes(v, {'double'}, {'scalar', 'integer', 'nonnegative'}); % degree v
validateattributes(u, {'double'}, {'scalar', 'integer'}); % order u
validateattributes(N, {'double'}, {'scalar', 'integer', 'nonnegative'}); % truncation order N
validateattributes(c, {'double'}, {'scalar'});
validateattributes(time_vec, {'double'}, {'row'}); % must be row vector
validateattributes(tau, {'double'}, {'scalar', 'nonnegative'});

%% Constant scaling factor
const = c/(4*pi*r*r_s); 

%% Find the time range during which intersection occurs
range = (time_vec >= abs(r-r_s)/c + tau) & (time_vec <= abs(r+r_s)/c + tau);

%% Find the values of Pvu(cos_theta_dash) 
% cos_theta_dash is the cosine of the elevation angle of the ring which is
% the intersection of the observation sphere and the expanding spherical
% surface measured w.r.t. the origin of the source
cos_theta_dash = find_ring_elev_wrt_source_delayed(r, r_s, c, time_vec, tau);

% Change values in the vicinity of +/-1
cos_theta_dash(cos_theta_dash < -1) = -1 + 1e-10;
cos_theta_dash(cos_theta_dash > 1) = 1 - 1e-10;

Pvu_mat = zeros(1, numel(time_vec)); % shape dimension
norm_Pvu_all = legendre(v, cos_theta_dash, 'norm'); % normalised associated Legendre polynomials

% extract normalised associated Legendre polynomial of order abs(u)
norm_Pvu = norm_Pvu_all(abs(u)+1, :); 

% modify Pvu according to the sign of order u
if u < 0
    % do nothing
else % if u >= 0
    norm_Pvu = (-1)^u * norm_Pvu;
end

Pvu_mat(range) = norm_Pvu(range);
%% Find the values of Pnu(cos_theta)
% cos_theta is the cosine of the elevation angle of the ring which is the
% intersection of the observation sphere and the expanding spherical
% surface measured w.r.t. the origin of the observation sphere
cos_theta = find_ring_elev_wrt_obv_delayed(r, r_s, c, time_vec, tau); 

cos_theta(cos_theta > 1) = 1-1e-10;
cos_theta(cos_theta < -1) = -1 + 1e-10;

Pnu_mat = zeros((N+1)^2, numel(time_vec)); % shape dimension

for n = abs(u):N % n must be greater than or equal to abs(u)
    norm_Pnu_all = legendre(n, cos_theta, 'norm'); % normalised associated Legendre polynomials
    
    % extract normalised associated Legendre polynomial of order abs(u)
    norm_Pnu = norm_Pnu_all(abs(u)+1, :); 
    
    % modify the value according to the sign of u
    if u < 0
        % do nothing
    else % if u >= 0
        norm_Pnu = (-1)^u * norm_Pnu;
    end
    
    Pnu_mat(n^2+n+1+u, range) = norm_Pnu(range); % allocate the correct row depending on the degree n
end
%% Find SH_coeff
SH_coeff = const*Pvu_mat.*Pnu_mat;
end