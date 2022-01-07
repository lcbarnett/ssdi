function [A,C,K,V] = transform_ss(A,C,K,V)

% Transform VAR model parameters to decorrelate and normalise residuals;
% i.e., so that residuals covariance matrix is the identity matrix.

n = size(C,1);
SQRTV  = chol(V,'lower');
ISQRTV = inv(SQRTV);
% A is unchanged
C = ISQRTV*C;
K = K*SQRTV;
V = eye(n);
