
function [n_mat, m_mat, nm_mask] = getDegreeOrderPairs(N)
% [n_mat, m_mat, nm_mask] = getDegreeOrderPairs(N)
%
% This function returns the degree-order pairs i.e. n-m pairs stored in 
% n_mat and m_mat for given spherical harmonic truncations order N. 
% The mask nm_mask masks out the unused n and m.
%
% Input:
% N - spherical harmonic truncation orders, must be a row vector
%
% Outputs:
% n_mat - degrees
% m_mat - orders
% nm_mask - mask to mask out unused n and m
% size(n_mat) = size(m_mat) = size(nm_mask) = [(max(N)+1)^2, size(N, 2)]

%% Check all elements in row vector N are nonnegative integers
validateattributes(N, {'double'}, {'integer', 'row', 'nonnegative'});

%% Calculate n and m
N_max = max(N);
n = repelem((0:N_max), 2.*(0:N_max)+1).'; % (N_max + 1)^2 by 1 column vector
m = (1:(N_max+1)^2).' - n.^2 - n - 1; % (N_max + 1)^2 by 1 column vector

n_mat = repmat(n, 1, size(N, 2));
m_mat = repmat(m, 1, size(N, 2));

double_nm_mask = zeros(size(n_mat));
for idx = 1:size(N, 2)
    row_limit = (N(idx)+1)^2;
    double_nm_mask(1:row_limit, idx) = 1;
end

nm_mask = logical(double_nm_mask);

n_mat = n_mat .* double_nm_mask;
m_mat = m_mat .* double_nm_mask;
end