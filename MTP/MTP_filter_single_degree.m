
function  y = MTP_filter_single_degree(n, r1, r2, c, t)
% y = MTP_filter_single_degree(n, r1, r2, c, t)
%
% This function calculates the MTP filter of a single degree n
% The MTP filter is the inverse Fourier transform of
% -i*k*jn(k*r1)*hn^(2)(k*r2)
%
% Inputs:
% n - degree of the MTP filter, must be a scalar
% r1 - small radius, in metres, must be a scalar
% r2 - big radius, in metres, must be a scalar
% c - speed of sound, in metres per second, scalar
% t - sampled time instances, in seconds, must be a row vector
%
% Output:
% y - The MTP filter of degree n evaluated at time samples in t
%     size(y) = [1, numel(t)]

%% Validate inputs
validateattributes(n, {'double'}, {'scalar', 'integer', 'nonnegative'});
validateattributes(c, {'double'}, {'scalar'});
validateattributes(t, {'double'}, {'row'});

validateattributes(r1, {'double'}, {'scalar', 'nonnegative'});
validateattributes(r2, {'double'}, {'scalar', 'nonnegative'});

% r2 must be greater than r1
if r2 <= r1
    error('@@ MTP_filter_single_degree: r2 must be greater than r1');
else
    % do nothing
end

%% Calculate the argument of Pn
y = zeros(size(t)); % shaping the dimension of output
cosine = zeros(size(t)); % argument of the Legendre function, dimension shaping

% Restrict the argument of the Legendre function to the time interval 
% abs(r1-r2) <= ct <= (r1+r2)
range = (t >= abs(r1-r2)/c) & (t <= abs(r1+r2)/c); 

% Argument of the Legendre function
cosine(range)= (r1^2 + r2^2 - c^2 * t(range).^2)/ (2*r1*r2);

cosine(cosine > 1-1e-10) = 1-1e-10;
cosine(cosine < -1+1e-10) = -1+1e-10;

% Scaling factor
scaling_factor = c/(2*r1*r2);

%% Calculate the Legendre function of degree n
% Use legendre instead of legendreP
assoc_Legendre = legendre(n, cosine); 
Legendre_func = assoc_Legendre(1, :)* scaling_factor;

y(range) = Legendre_func(range);
end
