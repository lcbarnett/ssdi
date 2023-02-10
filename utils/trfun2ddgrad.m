function [G,mG] = trfun2ddgrad(L,H)

% Calculate gradient of spectral dynamical dependence
% of projection L from transfer function H.
%
% NOTE 1: assumes identity residuals covariance matrix
% NOTE 2: projection L MUST be orthonormal!!!

h = size(H,3);
[n,m] = size(L);
g = zeros(n,m,h);
for k = 1:h % over [0,pi]
	Hk  = H(:,:,k);
	HLk = Hk'*L;
	g(:,:,k) = real((Hk*HLk)/(HLk'*HLk)); % grad/2
end

% Integrate frequency-domain derivative (trapezoidal rule) and
% subtract 2*L to get time-domain gradient (see note below)

G = sum(g(:,:,1:end-1)+g(:,:,2:end),3)/(h-1) - 2*L;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: the gradients are
%
% G <--- G-L*L'*G;   Grassmannian: Edelman-Arias-Smith, eq. (2.70)
% G <--- G-L*G'*L;   Stiefel:      Edelman-Arias-Smith, eq. (2.53)
%
% Here these turn out to be the same! In fact in our case
%
% L*L'*G = L*G'*L = 2*L
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&

if nargout > 1
	mG = sqrt(sum(G(:).^2)); % magnitude of gradient
end
