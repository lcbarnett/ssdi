function Q = orthonormalise(P)

% Returns an orthonormal basis Q for the range of P.
% That is, Q'*Q = I, the columns of Q span the same space as
% the columns of P, and the number of columns of Q is the
% rank of P.
%
% This method uses Singular Value Decompoistion; there are other
% (possibly more efficient) ways of doing this, such as QR decomposition,
% but SVD is very stable.

[Q,~] = svd(P,'econ');
