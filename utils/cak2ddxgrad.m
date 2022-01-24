function [G,mG] = cak2ddxgrad(L,CAK)

% Calculate gradient of proxy dynamical dependence
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

r = size(CAK,3);
n = size(L,1);
P = L*L';
g = zeros(n);
for k = 1:r
	Q = CAK(:,:,k);
	QT = Q';
	g = g + Q*QT - QT*P*Q - Q*P*QT;
end
G = 2*g*L;
G = G - P*G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: the gradients are
%
% G <--- G-L*L'*G;   Grassmannian: Edelman-Arias-Smith, eq. (2.70)
% G <--- G-L*G'*L;   Stiefel:      Edelman-Arias-Smith, eq. (2.53)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&

if nargout > 1
	mG = sqrt(sum(G(:).^2)); % magnitude of gradient
end
