function d = trfun2sdd(L,H)

% Calculate spectral dynamical dependence of projection L
% from transfer function H.
%
% NOTE 1: assumes identity residuals covariance matrix
% NOTE 2: projection L MUST be orthonormal!!!

h = size(H,3);
d = zeros(h,1);
LTL = L*L';
for k = 1:h % over [0,pi]
	LHk = L'*H(:,:,k);
	LHLTk = LHk*L;
    d(k) = logdet(LHk*LHk') - logdet(LHLTk*LHLTk');
end
