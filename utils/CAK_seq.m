function U = CAK_seq(A,C,K)

% Return the sequence U_k = CA^{k-1}K, k = 1,...,r of n times matrices

r = size(A,1);
n = size(C,1);
Ak = eye(r);
U = zeros(n,n,r);
U(:,:,1) = C*K;
for k = 2:r
	Ak = Ak*A;
	U(:,:,k) = C*Ak*K;
end
