function L = Q2L(Q,m)

% Returns orthonormal bases L corresponding to involutions Q
%
% Could use QR decomposition, but SVD may be more stable.
%
% Ref: Z. Lai, L.-H. Lim and  K.Ye, "Simpler Grassmannian optimization"
% arXiv:2009.13502 [math.OC] - https://arxiv.org/abs/2009.13502

[n,~,r] = size(Q);
L = zeros(n,m,r);
for k = 1:r
	[U,~] = svd((eye(n)+Q(:,:,k))/2);
%	[U,~] = qr((eye(n)+Q(:,:,k))/2); % computationally cheaper than SVD
	L(:,:,k) = U(:,1:m);
end
