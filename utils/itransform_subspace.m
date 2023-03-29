function L0 = itransform_subspace(L,V0)

% Inverse-transform projections (subspaces) back for correlated residuals
%
% L        transformed orthonormal subspace basis
% V0       untransformed residuals covariance matrix
% L0       untransformed orthonormal subspace basis

assert(ndims(L) >= 2);
siz = size(L);
xdims = prod(siz(3:end)); % extra dimesions
IV0LCT = inv(chol(V0,'lower')');
L0 = reshape(L,siz(1),siz(2),xdims);
for k = 1:xdims
	L0(:,:,k) = orthonormalise(IV0LCT*L0(:,:,k)); % note that IV0LCT'*L will NOT generally be orthonormal!
end
L0 = reshape(L0,siz);
