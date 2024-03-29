function [D,d] = trfun2dd(L,H)

% Calculate spectral dynamical dependence of projection L
% from transfer function H.
%
% L    orthonormal subspace basis
% H    transfer function
%
% NOTE: assumes identity residuals covariance matrix

h = size(H,3);
d = zeros(h,1);
for k = 1:h % over [0,pi]
	Qk = H(:,:,k)'*L;
	d(k) = sum(log(diag(chol(Qk'*Qk)))); % (log-determinant)/2
end
D = sum(d(1:end-1)+d(2:end))/(h-1); % integrate frequency-domain DD (trapezoidal rule) to get time-domain DD
