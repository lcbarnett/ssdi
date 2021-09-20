function D = sliss2dd(L,A,C,K,ldwork)

% Calculate dynamical dependence of projection L for
% innovations-form state-space model with parameters A,C,K.
%
% Like iss2dd, but uses SLICOT DARE solver (slidare from MVGC2)
% See opt_es_sldd for 'ldwork'.
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L MUST be orthonormal!!!

% Calculate residuals covariance matrix V of projected model (solve DARE)

[~,info,V] = slidare(A',C'*L,K*K',[],(K*L)',ldwork);

if info.rep ~= 0 % DARE failed
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
