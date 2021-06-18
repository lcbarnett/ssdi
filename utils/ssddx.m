function D = ssddx(L,M,CAK)

% Calculate proxy dynamical dependence of projection L for
% innovations-form state-space model with parameters A,C,K,
% for preoptimisation
%
% U_k, k = 1,...,r is the sequence of m x (n-m) matrices CA^{k-1}K
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L, M MUST be orthonormal and orthogonal!!!

% Calculate Frobenius norm of dynamical independence condition

r = size(CAK,3);
D = 0;
for k = 1:r
	Wk = L'*CAK(:,:,k)*M;
	Wk2 = Wk.*Wk;
	D = D + sum(Wk2(:));
end
