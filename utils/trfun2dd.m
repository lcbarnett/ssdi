function [D,d] = trfun2dd(L,H)

% Calculate spectral dynamical dependence of projection L
% from transfer function H.
%
% NOTE 1: assumes identity residuals covariance matrix
% NOTE 2: projection L MUST be orthonormal!!!

h = size(H,3);
d = zeros(h,1);
LT = L';
for k = 1:h % over [0,pi]
	LHk = LT*H(:,:,k);
    d(k) = logdet(LHk*LHk');
end

D = trapz(d)/(h-1); % integrate frequency-domain DD to get time-domain DD
