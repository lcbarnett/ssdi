function CAK = iss2cak(A,C,K)

% Return the sequence CA^{k-1}K, k = 1,...,r of n x n matrices

r = size(A,1);
n = size(C,1);
Ak = eye(r);
CAK = zeros(n,n,r);
CAK(:,:,1) = C*K;
for k = 2:r
	Ak = Ak*A;
	CAK(:,:,k) = C*Ak*K;
end
