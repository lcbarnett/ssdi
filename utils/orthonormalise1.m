function [L,M] = orthonormalise1(P)

% Returns an orthonormal basis L for the range of P and a basis
% M for the orthoganal subspace. That is, L'*L = I, M'*M = I ,
% L'*M = 0, the columns of L span the same space as the columns
% of P, while the columns of L span the orthogonal subspace.

m = size(P,2);
[U,~] = svd(P);
L = U(:,1:m);
M = U(:,m+1:end);
