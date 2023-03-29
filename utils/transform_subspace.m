function L = transform_subspace(L0,V0)

% Transform projections (subspaces) for correlated residuals
%
% L0       untransformed orthonormal subspace basis
% V0       untransformed residuals covariance matrix
% L        transformed orthonormal subspace basis

assert(ndims(L0) >= 2);
siz = size(L0);
xdims = prod(siz(3:end)); % extra dimesions
V0LCT = chol(V0,'lower')';
L = reshape(L0,siz(1),siz(2),xdims);
for k = 1:xdims
	L(:,:,k) = orthonormalise(V0LCT*L(:,:,k)); % note that V0LCT'*L0 will NOT generally be orthonormal!
end
L = reshape(L,siz);
