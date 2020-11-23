function L = orthonormalise(P)

% Returns an orthonormal basis L for the range of P.
% That is, L'*L = I, the columns of L span the same space as
% the columns of P, and the number of columns of L is the
% rank of P.
%
% This method uses Singular Value Decompoistion; there are other
% (possibly more efficient) ways of doing this, such as QR decomposition,
% but SVD is very stable.

[L,~] = svd(P,'econ');
