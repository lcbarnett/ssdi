function [CE,DD,CI] = iss2ce(L,A,C,K,V,G,P)

% Calculate Causal Emergence (CE) of projection L for
% innovations-form state-space model with parameters A,C,K,V
%
% L        orthonormal subspace basis
% A,C,K,V  ISS parameters
% G        ISS covariance matrix (if empty, calculate)
% P        per-element DARE solutions (if empty, calculate)
%
% NOTE: DOES NOT assume identity residuals covariance matrix!

[n,m] = size(C);

VSR = chol(V,'lower');

if isempty(G) || isempty(P)
	KV = K*VSR;
	KVK = KV*KV';
end

% First term

if isempty(P)
	P = zeros(m,m,n);
	for i = 1:n
		[~,~,~,~,P(:,:,i)] = mdare(A,KVK,V(i,i),K*V(:,i))
	end
end
ce1 = zeros(n,1);
for i = 1:n
	CPi = C*chol(P(:,:,i),'lower');
	ce1(i) = sum(log(diag(L'chol(CPi*CPi'+V,'lower'))));
end
CE1 = 2*sum(ce1);

% Second term

LV = L'*VSR;
[~,~,~,~,PL] = mdare(A,KVK,LV*LV',K*V*L);
CPL = C*chol(PL,'lower');
CE2 = 2*sum(log(diag(L'*chol(CPL*CPL'+V))));

% Third term

if isempty(G)
	CM = C*chol(dlyap(A,KVK),'lower');
	G  = CM*CM' + V;
end
CE3 = (n-1)*2*sum(log(diag(L'*chol(G,'lower'))));

% Causal Emergence

CE = CE1 - CE2 - CE3;

if nargout > 1

	CE4 = sum(log(diag(LV)));

	% Dynamical Dependence

	DD = CE2 - CE4;

	if nargout > 2

		% Co-information

		CI = CE1 - CE4 - CE3;

	end
end
