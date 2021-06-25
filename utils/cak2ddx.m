function D = cak2ddx(L,CAK)

% Calculate proxy dynamical dependence of projection L for
% innovations-form state-space model with parameters A,C,K,
% for preoptimisation
%
% CAK_k, k = 1,...,r is the sequence of m x (n-m) matrices CA^{k-1}K
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L orthogonal!!!

% Calculate Frobenius norm of dynamical independence condition

r = size(CAK,3);
D = 0;
for k = 1:r
	LCAKk = L'*CAK(:,:,k);
	LCAKLTk = LCAKk*L;
	D1k = LCAKk.*LCAKk;
	D2k = LCAKLTk.*LCAKLTk;
	D = D + sum(D1k(:)) - sum(D2k(:));
end