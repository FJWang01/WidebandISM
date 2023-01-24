function Ylm=sphHarm(l,m,theta,phi)
% This function calculates the spherical harmonics function of degree l and order m at location (theta, phi)
% Input
% l - degree
% m - order
% theta, phi - column vectors that contain the elevation and azimuth angles of
% locations of interest, in radians
if ~isequal(size(theta),size(phi))
	error('@@ sphHarm: theta and phi need to be the same size')
end

Qlm=(-1)^m*sqrt((2*l+1)/(8*pi)); % (-1)^m Q_l 
Ylm=zeros(size(theta)); % defines shape and init for output
S_max=numel(theta);
for s=1:S_max % walk through theta,phi matrices as vectors
	Sl=legendre(l,cos(theta(s)),'sch'); % vector in m, scalar in theta
	Slm=Sl(abs(m)+1,1); % pull out scalar S_l^|m| (scalar theta)
	Ylm(s)=Qlm*Slm*exp(1j*m*phi(s)); % scalar (implicit 'conj' for m<0)
end

if m==0
	Ylm=sqrt(2)*Ylm; % m=0 adjustment
elseif m<0
	Ylm=(-1)^m*Ylm; % m<0 adjustment
end