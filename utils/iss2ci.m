function CI = iss2ci(L,A,C,K,V,G,P)

% Calculate co-information (for CE) of projection L for
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
ci1 = zeros(n,1);
for i = 1:n
	CPi = C*chol(P(:,:,i),'lower');
	ci1(i) = sum(log(diag(L'chol(CPi*CPi'+V,'lower'))));
end
CI1 = 2*sum(ci1);

% Second term

CI2 = 2*sum(log(diag(L'*VSR)));

% Third term

if isempty(G)
	CM = C*chol(dlyap(A,KVK),'lower');
	G  = CM*CM' + V;
end
CI3 = (n-1)*2*sum(log(diag(L'*chol(G,'lower'))));

% Co-information

CI = CI1 - CI2 - CI3;
