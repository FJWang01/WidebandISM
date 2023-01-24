function y_bes = sph_Bessel_2(n, x)
% This function calculates the spherical Bessel function of the 2nd kind
% of order n at argument x
% 
% Input
% n - order of the spherical Bessel function of the 2nd kind
%     n must be non-negative integer
% x - argument of the function, must be real
%
% Note:
% x and n must be of the same size
%
% Output
% y_bes - values of the spherical Bessel function of the 2nd kind of order 
%         n evaluated at x
%         size(y_bes) = size(n) = size(x)

%% Notice:
% The values of the spherical Bessel function of the second kind with 
% negative arguments are calculated using the relationship:
% yn(-x) = (-1)^(n+1)*yn(x) ---- (1)
% See equation 10.47.14 in https://dlmf.nist.gov/10.47

% The order n must be non-negative integer, see section 10.47(i)

% Also, see Wolfram Alpha
% https://mathworld.wolfram.com/SphericalBesselFunctionoftheSecondKind.html
% https://reference.wolfram.com/language/ref/SphericalBesselY.html

% Another example 
% https://math.stackexchange.com/questions/170013/value-of-a-scaled-bessel-function-for-negative-argument
%% Check if nu and x are of the same size
if ~isequal(size(n),size(x))
	error('@@ sph_Bessel_2: n and x must be of the same size');
else
    % do nothing
end
%% n must be non-negative integer
validateattributes(n, {'double'}, {'integer', 'nonnegative'});
%% Calculate yn(abs(x))
% Formula of yn(x)
% yn(x) = sqrt(pi./(2x)).*bessely(n+0.5, x)
Y = bessely(n+0.5,abs(x)); % Bessel function of second kind with order n+0.5
y_bes = sqrt(pi./(2*abs(x))).*Y;
%% Use (1) to find values of yn(x) when x is negative
y_bes(x<0) = (-1).^(n(x<0)+1).*y_bes(x<0); % apply (1)
end