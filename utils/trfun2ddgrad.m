function G = trfun2ddgrad(L,H)

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
	Qk = Hk'*L;
	g(:,:,k) = real((Hk*Qk)/(Qk'*Qk)); % grad/2
end
G = sum(g(:,:,1:end-1)+g(:,:,2:end),3)/(h-1); % integrate frequency-domain gradient (trapezoidal rule) to get time-domain gradient
G = G-(L*L')*G; % Edelman-Arias-Smith, eq. (2.70)
