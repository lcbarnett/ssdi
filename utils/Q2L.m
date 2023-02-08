function L = Q2L(Q,m)

% Returns orthonormal bases L corresponding to involutions Q
%
% Ref: Z. Lai, L.-H. Lim and  K.Ye, "Simpler Grassmannian optimization"
% arXiv:2009.13502 [math.OC] - https://arxiv.org/abs/2009.13502

[n,~,r] = size(Q);
L = zeros(n,m,r);
for k = 1:r
	[U,~] = svd((eye(n)+Q(:,:,k))/2);
	L(:,:,k) = U(:,1:m);
end
