function [A,V] = transform_var(A,V)

% Transform VAR model parameters to decorrelate and normalise residuals;
% i.e., so that residuals covariance matrix is the identity matrix.

[n,~,p] = size(A);
SQRTV  = chol(V,'lower');
ISQRTV = inv(SQRTV);
for k = 1:p
	A(:,:,k) = ISQRTV*A(:,:,k)*SQRTV;
end
V = eye(n);
