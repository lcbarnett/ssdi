function D = cak2ddx(L,CAK)

% Calculate proxy dynamical dependence
%
% For an innovations-form state-space model with parameters (A,C,K)
%
%     CAK_k, k = 1,...,r is the sequence of n x n matrices CA^{k-1}K
%
% where r is the ISS model order; see iss2cak.m.
%
% For a VAR model with coefficients sequence A, we may supply simply CAK = A.
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L orthogonal!!!

% Calculate Frobenius norm of dynamical independence condition

r = size(CAK,3);
D = 0;
for k = 1:r
	LCAKk = L'*CAK(:,:,k);
	LCAKLTk = LCAKk*L;
	D1k = LCAKk.^2;
	D2k = LCAKLTk.^2;
	D = D + sum(D1k(:)) - sum(D2k(:));
end
