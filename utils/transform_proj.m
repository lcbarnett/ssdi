function L = transform_proj(L,V)

% Inverse-transform projections (hyperplanes) for correlated residuals

assert(ndims(L) >= 2);
siz = size(L);
xdims = prod(siz(3:end)); % extra dimesions
L = reshape(L,siz(1),siz(2),xdims);
SQRTV  = chol(V,'lower');
ISQRTV = inv(SQRTV);
for k = 1:xdims
	L(:,:,k) = ISQRTV'*L(:,:,k);
end
L = reshape(L,siz);
