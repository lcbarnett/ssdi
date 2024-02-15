function [G,P] = iss2ce_precomp(A,C,K,V)

% Precompute observed process covariance matrix G and per-element state prediction
% error covariance matrices P for Causal Emergence computation; see iss2ce.m.
%
% A,C,K,V  ISS parameters
%
% NOTE: does NOT assume normalisation of residuals covariance matrix V.

[n,m] = size(C);

VCHOL = chol(V,'lower');
KV  = K*VCHOL;
KVK = KV*KV';

M = dlyap(A,KVK); % state covariance
G = C*dlyap(A,KVK)*C' + V;
P = zeros(m,m,n);
for i = 1:n
	[~,~,~,~,P(:,:,i)] = mdare(A,C(i,:),KVK,V(i,i),K*V(:,i));
end
