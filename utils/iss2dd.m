function D = iss2dd(L,A,C,K)

% Calculate dynamical dependence of projection L for
% innovations-form state-space model with parameters A,C,K.
%
% L        orthonormal subspace basis
% A,C,K    ISS parameters
%
% NOTE: assumes identity residuals covariance matrix

% Calculate residuals covariance matrix V of projected model (solve DARE)

[~,V,rep] = mdare(A,L'*C,K*K',[],K*L); % mdare from MVGC2

if rep < 0  || rep > 1e-08 % DARE failed
	D = NaN;
	return
end

% D = log-determinant of residuals covariance matrix V

[R,p] = chol(V);
if p == 0
	D = 2*sum(log(diag(R)));
else
	D = NaN; % fail: V not positive-definite
end
