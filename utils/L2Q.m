function Q = L2Q(L)

% Returns involutions Q corresponding to orthonormal bases L
%
% Ref: Z. Lai, L.-H. Lim and  K.Ye, "Simpler Grassmannian optimization"
% arXiv:2009.13502 [math.OC] - https://arxiv.org/abs/2009.13502

[n,m,r] = size(L);
Q = zeros(n,n,r);
for k = 1:r
	Q(:,:,k) = 2*L*L'-eye(n);
end
