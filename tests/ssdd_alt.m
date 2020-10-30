function [D,E] = ssdd_alt(P,Q,A,C,K)

% Calculate dynamical dependence of projection P for
% innovations-form state-space model with parameters A,C,K.
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection P MUST be orthonormal!!!

% Calculate residuals covariance matrix V of projected model (solve DARE)

[~,V,rep] = ssdare(A,P'*C,K*K',eye(size(P,2)),K*P);

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

r = size(A,1);
[n,m] = size(P);
PC = P'*C;
KQ = K*Q;
%A = A-K*C;
%A = A/max(abs(eig(A)));
J = zeros(m*m,n-m);
Ak = eye(r);
J(1:m,:) = PC*Ak*KQ;
for k = 1:r-1
	Ak = A*Ak;
	J(k*m+1:(k+1)*m,:) = PC*Ak*KQ;
end
E = mean(J(:).^2);
%E = mean([PC(:);KQ(:)].^2);

%{
r = size(A,1);
PC = P'*C;
KQ = K*Q;
An = A/max(abs(eig(A)));
Ak = eye(r);
E = norm(PC*Ak*KQ);
for k = 1:r-1
	Ak = An*Ak;
	E = E + norm(PC*Ak*KQ);
end
E = E/r;
%}

%{
r = size(A,1);
PC = P'*C;
KQ = K*Q;
B = A-K*C;
Bn = B/max(abs(eig(B)));
Bk = eye(r);
E = norm(PC*Bk*KQ);
for k = 1:r-1
	Bk = Bn*Bk;
	E = E + norm(PC*Bk*KQ);
end
E = E/r;
%}
