
function j_bes = sph_Bessel_1(n, x)
% This function calculates the spherical Bessel function of the 1st kind
% of order n with argument x
% 
% Input
% n - order of the spherical Bessel function of the 1st kind
%     n must be non-negative integer
% x - argument of the function, must be real
% 
% x and n must be of the same size
%
% Output
% j_bes - values of the spherical Bessel function of the 1st kind of order 
%         n evaluated at x
% size(j_bes) = size(n)
%
% Notice:
% The values of the spherical Bessel function of the first kind with 
% negative arguments are calculated using the relationship:
% jn(-x) = (-1)^n*jn(x) ---- (1)
% See equation 10.47.14 in https://dlmf.nist.gov/10.47
% The order n must be a non-negative integer, see section 10.47(i)
%
% Also, see Wolfram Alpha
% https://mathworld.wolfram.com/SphericalBesselFunctionoftheFirstKind.html
% https://reference.wolfram.com/language/ref/SphericalBesselJ.html
%
% Another explanation 
% https://math.stackexchange.com/questions/170013/value-of-a-scaled-bessel-function-for-negative-argument

%% Check if n and x are of the same size
if ~isequal(size(n),size(x))
	error('@@ sph_Bessel_1: n and x must be of the same size');
else
    % do nothing
end

%% Check if any n is negative
validateattributes(n, {'double'}, {'integer', 'nonnegative'});

%% Calculate jn(abs(x))
% Formula:
% jn(x) = sqrt(pi./(2x)).*besselj(n+0.5, x), x must be positive
J = besselj(n+0.5, abs(x)); % Bessel function of first kind of order (n+0.5)
j_bes = sqrt(pi./(2*abs(x))).*J; % jn(abs(x))

%% Use equation (1) to find the values of jn(x) when x<0
j_bes(x<0) = (-1).^n(x<0).*j_bes(x<0); % apply (1)

%% Fix the values at x == 0
L_n0_x0 = (n == 0) & (x == 0); % L_n0_x0 = 1 when n == 0 and x == 0
j_bes(L_n0_x0) = 1; % j_0(0) = 1

L_x0 = (n > 0) & (x == 0); % L_x0 = 1 when n ~= 0 and x == 0
j_bes(L_x0) = 0; % jn(0) = 0 when n ~= 0
end
