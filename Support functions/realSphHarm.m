function y = realSphHarm(n, m, theta, phi)
% This function calculates the real spherical harmonic function of degree n
% and order m at elevation angle(s) theta and azimuth angle(s) phi
% Input
% n - degree, scalar
% m - order, scalar
% theta - elevation angle, in radians
% phi - azimuth angle, in radians

% Output
% y - the real spherical harmonic function of degree n and order m 
% at elevation angle theta and azimuth angle phi

%% theta and phi must be of same size
if ~isequal(size(theta),size(phi))
	error('@@ realSphHarm: theta and phi need to be the same size');
end

%% perform calculations
y = zeros(size(theta)); % shape the dimension of output
S_max = numel(theta);

for s = 1:S_max % go through each theta-phi pair
    Pnm_all = legendre(n, cos(theta(s)), 'norm');
    Pnm = Pnm_all(abs(m)+1, 1);
    if m < 0
        K = sin(abs(m)*phi(s));
    elseif m > 0
        K = cos(m*phi(s));
    else % m = 0
        K = 1/sqrt(2);
    end 
    y(s) = 1/sqrt(pi)*K*Pnm;
end
end