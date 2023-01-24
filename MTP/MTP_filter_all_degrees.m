
function  y_all_degrees = MTP_filter_all_degrees(n, r1, r2, c, t)
% y_all_degrees = MTP_filter_all_degrees(n, r1, r2, c, t)
%
% This function calculates the MTP filters of all degrees in n
% The MTP filter is the inverse Fourier transform of
% -i*k*jn(k*r1)*hn^(2)(k*r2) 
%
% Inputs:
% n - degrees of MTP filters, must be a column vector
% r1 - small radius, in metres, must be a scalar
% r2 - big radius, in metres, must be a scalar
% c - speed of sound, in metres per second, scalar
% t - sampled time instances, in seconds, must be a row vector
%
% Output:
% y_all_degrees - MTP filters of all degrees in n evaluated at time samples 
%                 in t
%                 size(y_all_degrees) = [numel(n), numel(t)]
%

%% Validate inputs
validateattributes(n, {'double'}, {'column', 'integer', 'nonnegative'});
validateattributes(r1, {'double'}, {'scalar', 'nonnegative'});
validateattributes(r2, {'double'}, {'scalar', 'nonnegative'});
validateattributes(c, {'double'}, {'scalar'});
validateattributes(t, {'double'}, {'row'});

% r2 must be greater than r1
if r2 <= r1
    error('@@ MTP_filter_all_degrees: r2 must be greater than r1');
else
    % do nothing
end

%% Calculate the argument of Pn
cosine = zeros(1, numel(t)); % argument of the Legendre function, dimension shaping

% Restrict the argument of the Legendre function to the time interval 
% abs(r1-r2) <= ct <= (r1+r2)
range = (t >= abs(r1-r2)/c) & (t <= abs(r1+r2)/c); 

% Argument of the Legendre function
cosine(range)= (r1^2 + r2^2 - c^2 * t(range).^2)/ (2*r1*r2);

cosine(cosine > 1-1e-10) = 1-1e-10;
cosine(cosine < -1+1e-10) = -1+1e-10;

% Scaling factor
scaling_factor = c/(2*r1*r2);

%% Calculate the Legendre functions of degrees in n
% Use legendre instead of legendreP
% Find unique values in n
[unique_n, ~, ic] = unique(n, 'stable');

% Calculate Legendre function for each unique n value
y = zeros(numel(unique_n), numel(t)); % shaping the dimension of output

% Go through each unique n value
for idx = 1:numel(unique_n)
    assoc_Legendre = legendre(unique_n(idx), cosine); 
    Legendre_func = assoc_Legendre(1, :)* scaling_factor;
    y(idx, range) = Legendre_func(range);
end

y_all_degrees = y(ic, :);
end
