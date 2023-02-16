function [L,M] = orthonormalise(X)

% Returns an orthonormal basis L for the range of X.
% That is, L'*L = I, the columns of L span the same space as
% the columns of X, and the number of columns of L is the
% rank of X. Optionally, also return orthonormal basis M for
% orthogonal subspace.
%
% Could use QR decomposition, but SVD may be more stable.

if nargout < 2
	[L,~] = svd(X,0);
%	[L,~] = qr(X,0);
else
	[m,n] = size(X);
	[U,~] = svd(X);
%	[U,~] = qr(X);
	L = U(:,1:m);
	M = U(:,m+1:n);
end
