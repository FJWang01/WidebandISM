function DM=WignerDM(N,alpha,beta,gamma)
% build Wigner D matrix (eq. (20), Ben Hagai Nov 2012) 
% first, gamma about z, then beta about y, then alpha about z again (counterclocwise)
DM=zeros((N+1)^2);
    for n=0:N
        for m=-n:n
            DM(n^2+n+m+1,n^2+1:(n+1)^2)=...
                wignerd(n,m,alpha,beta,gamma);
        end
    end
end

function D=wignerd(n,m,alpha,beta,gamma)
% wignerd.m
% ------------
%       Get Wigner-D Coefficients.
%       Based on Rafaely 2008, equations (11-12).
%
% Syntax
% ------------
%     D=wignerd(n,m,m2,alpha,beta,gamma)
%
% Input
% ------------
%     Required
%           n,m - scalar - the spherical indices
%           alpha,beta,gamma - Euler angles
%
% Output
% ------------
%         D - defined as D_{mm'}^n, where m'=-n:n
%
% Created/Modified by
% ------------
%     Ilan Ben Hagai, 1-Nov-2010

%%
    m2=-n:n;
    
    epsilon=1.*(m2>=m) + (-1).^(m2-m).*(m2<m);
    mu=abs(m-m2);
    nu=abs(m+m2);
    s=n-(mu+nu)/2;
    Ps=zeros(1,numel(mu));
    for mIdx=1:numel(mu)
        curMu=mu(mIdx);
        curNu=nu(mIdx);
        curS=s(mIdx);
        polyCoeffs=orth_poly('Jacobi',curS,curMu,curNu);
        Ps(mIdx)= polyval(polyCoeffs,cos(beta));
    end
    
    % calculate the Wigner-d function (eq.12)
    d=epsilon.*sqrt(factorial(s).*factorial(s+mu+nu)./(factorial(s+mu).*factorial(s+nu))).*sin(beta/2).^mu.*cos(beta/2).^nu.*Ps;

    % Calculate the coefficients (eq.11) :
    D=exp(-1i*m*alpha-1i*m2*gamma).*d;
end

function pn=orth_poly(class,n,alpha,beta)
% generates an orthogonal polynomial

    if (nargin<4)||isempty(beta)
        beta=0;
    end
    if (nargin<3)||isempty(alpha);
        alpha=0;
    end

    % initialize (-1)'th and zero'th order polynomials
    pn=[];
    pnp1=1;

    for i=0:n
        pnm1=pn;
        pn=pnp1;
        switch class
            case 'Legendre'
                pnp1=((2*i+1)*[pn,0] - i*[0,0,pnm1])/(i+1);
            case 'Hermite'
                pnp1=2*[pn,0] - 2*i*[0,0,pnm1];
            case 'Laguerre'
                pnp1=((2*i+alpha+1)*[0,pn] -[pn,0] - (i+alpha)*[0,0,pnm1])/(i+1);
            case 'Jacobi'
                if (alpha~=0)||(beta~=0)
                    a1n=2*(i+1)*(i+alpha+beta+1)*(2*i+alpha+beta);
                    a2n=(2*i+alpha+beta+1)*(alpha^2-beta^2);
                    if (2*i+alpha+beta)<=150
                        a3n=gamma(2*i+alpha+beta+3)./gamma(2*i+alpha+beta);
                    else
                        a3n=exp(gammaln(2*i+alpha+beta+3)-gammaln(2*i+alpha+beta));
                    end
                    a4n=2*(i+alpha)*(i+beta)*(2*i+alpha+beta+2);
                    pnp1=(a2n*[0,pn] + a3n*[pn,0] - a4n*[0,0,pnm1])./a1n;
                else
                    pnp1=((2*i+1)*[pn,0] - i*[0,0,pnm1])/(i+1);
                end
        end

    end
end
