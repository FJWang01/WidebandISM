function h_Han_2 = sph_Hankel_2(n, x)
% This function calculates the spherical Hankel function of the second kind
% of order n at argument x
% 
% Input
% n - order of the spherical Hankel function of the second kind
%     n must be non-negative integer
% x - argument of the function, must be real
%   
% Note:
% n and x must be of the same size
%
% Output
% h_Han_2 - values of the spherical Hankel function of the second kind of 
%           order n evaluated at x
%           size(h_Han_2) = size(n) = size(x)

%% Check if n and x are of the same size
if ~isequal(size(n),size(x))
	error('@@ sph_Hankel_2: n and x must be of the same size');
else
    % do nothing
end
%% n must be non-negative integer
validateattributes(n, {'double'}, {'integer', 'nonnegative'});
%% Perform calculation using formula
% hn^(2)(x) = jn(x) - i*yn(x)
j_bes = sph_Bessel_1(n, x);
y_bes = sph_Bessel_2(n, x);
h_Han_2 = j_bes - 1i*y_bes; 
end