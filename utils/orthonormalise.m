function [L,M] = orthonormalise(P)

% Returns an orthonormal basis L for the range of P.
% That is, L'*L = I, the columns of L span the same space as
% the columns of P, and the number of columns of L is the
% rank of P. Optionally, also return orthonormal basis M for
% orthogonal subspace.
%
% This method uses Singular Value Decompoistion; there are other
% (possibly more efficient) ways of doing this, such as QR decomposition,
% but SVD is very stable.

if nargout < 2
	[L,~] = svd(P,'econ');
else
	m = size(P,2);
	[U,~] = svd(P);
	L = U(:,1:m);
	M = U(:,m+1:end);
end
