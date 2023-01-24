function h_Han_1 = sph_Hankel_1(nu, x)
% This function calculates the spherical Hankel function of the first kind
% with order nu and argument x
% I/O
% Input
% nu - order of the spherical Hankel function of the first kind
%    - elements of nu MUST be non-negative integer
% x - argument of the function, must be real

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warning: nu and x MUST be of the same size 
% 
% If one of them is scalar, please kindly use
% meshgrid to convert them to same size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output
% h_Han_1 - values of the spherical Hankel function of the first kind with order nu
% evaluated at x

%% Notice:
% The values of the spherical Hankel function of the first kind with 
% negative arguments satisfy the relationship:
% hn^(1)(-x) = (-1)^n*hn^(2)(x)
% See equation 10.47.15 in https://dlmf.nist.gov/10.47

% The order n MUST be a non-negative integer, see section 10.47(i)

% Also see Wolfram Alpha
% https://mathworld.wolfram.com/SphericalHankelFunctionoftheFirstKind.html
% https://reference.wolfram.com/language/ref/SphericalHankelH1.html
% Another explanation see 
% https://math.stackexchange.com/questions/170013/value-of-a-scaled-bessel-function-for-negative-argument
%% Check if nu and x are of the same size
if ~isequal(size(nu),size(x))
	error('@@ sph_Hankel_1: nu and x MUST be the same size')
else
    % do nothing
end
%% Check if any nu is negative
if ~isempty(find(nu<0,1))
    error('@@ sph_Hankel_1: elements of nu MUST be positive integer')
else
    % do nothing
end
%% Perform calculation using formula
% hn^(1)(x) = jn(x) + i*yn(x)
j_bes = sph_Bessel_1(nu, x);
y_bes = sph_Bessel_2(nu, x);
h_Han_1 = j_bes + 1i*y_bes; 
end