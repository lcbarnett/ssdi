function mad = dds_check(A,C,K,H,m,nsics)

% Check accuracy of spectral DD against state-space DD
% at given scale m across a set of random projections.

n = size(C,1);
L = rand_orthonormal(n,m,nsics); % (orthonormalised) random linear projections
D1 = zeros(nsics,1);
for k = 1:nsics
	D1(k) = iss2dd(L(:,:,k),A,C,K);
end
D2 = zeros(nsics,1);
for k = 1:nsics
	D2(k) = trfun2dd(L(:,:,k),H);
end
mad = maxabs(D1-D2);
