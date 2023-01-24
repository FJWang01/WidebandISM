
function PW_y_all_degrees = PW_MTP_filter_all_degrees(n, r, c, t)
% PW_y_all_degrees = PW_MTP_filter_all_degrees(n, r, c, t)
%
% This function calculates the plane wave MTP filters of degrees in n
% The plane wave MTP filter is the inverse Fourier transform of
% 4*pi*i^n*jn(kr) 
%
% Inputs:
% n - degrees of plane wave MTP filters, must be a column vector
% r - radius, in metres, must be a scalar
% c - speed of sound, in metres per second, scalar
% t - time samples, in seconds, must be a row vector
%
% Output:
% PW_y_all_degrees - plane wave MTP filters of degrees in n at sample times
%                    in t
%                    size(PW_y_all_degrees) = [numel(n), numel(t)]
%

%% Validate inputs
validateattributes(n, {'double'}, {'column', 'integer', 'nonnegative'});
validateattributes(r, {'double'}, {'scalar', 'nonnegative'});
validateattributes(c, {'double'}, {'scalar'});
validateattributes(t, {'double'}, {'row'});

%% Calculate the argument of the Legendre function
cosine = zeros(1, numel(t)); % argument of the Legendre function, dimension shaping

% Restrict the argument of the Legendre function to the time interval 
% -r/c <= t <= r/c
range = (t >= -r/c) & (t <= r/c); 

% Argument of the Legendre function
cosine(range)= - c * t(range)/ r;

cosine(cosine > 1-1e-10) = 1-1e-10;
cosine(cosine < -1+1e-10) = -1+1e-10;

% Scaling factor
scaling_factor = 2*pi*c/r;

%% Calculate the Legendre function of degree n
% use legendre instead of legendreP
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

PW_y_all_degrees = y(ic, :);
end
