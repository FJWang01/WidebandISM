
function Pnm = assoc_Legendre_poly(n, m, t)
% This function calculates the associated Legendre polynomial of degree n
% and order m at argument t
% Input
% n - degree, must be non-negative integer
% m - order, must be integer
% t - argument, between -1 and 1. Must be a row vector
%
% n and m must be of the same size
%
% Output
% Pnm - the associated Legendre polynomial of degree n
%       and order m at argument t
% size(Pnm) = [size(n, 1), size (n, 2), numel(t)]

%% Check if n and m are of the same size
if ~isequal(size(n),size(m))
	error('@@ assoc_Legendre_poly: n and m must be of the same size')
else
    % do nothing
end

%% Check if t is a row vector
validateattributes(t, {'double'}, {'nrows', 1}); % t must be a row vector
%% Mains
% Convert inputs into column vectors
n_col = n(:);
m_col = m(:);

Pnm_2D = zeros(numel(n_col), numel(t));
for pair_idx = 1:numel(n_col) % go through each degree order pair
    n_temp = n_col(pair_idx);
    m_temp = m_col(pair_idx);
    Pnm_all = legendre(n_temp, t);
    Pnm_desired = Pnm_all(abs(m_temp)+1, :);
    if m_temp < 0 % m<0 adjustment
        Pnm_desired = (-1)^abs(m_temp) * factorial(n_temp - abs(m_temp)) / factorial(n_temp + abs(m_temp)) * Pnm_desired;
    else
        % do nothing
    end
    Pnm_2D(pair_idx, :) = Pnm_desired;
end

% Rearrange Pnm_2D so that each page of Pnm has dimension == size(n)
Pnm = zeros(size(n, 1), size(n, 2), numel(t));
for arg_idx = 1:numel(t)
    Pnm(:, :, arg_idx) = reshape(Pnm_2D(:, arg_idx), size(n));
end
end