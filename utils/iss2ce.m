function [CI,DD] = iss2ce(L,A,C,K,V,G,P)

% Calculate co-information (CI) and Dynamical Dependence (DD) of projection L for
% innovations-form state-space model with parameters A,C,K,V. Causal Emergence (CE)
% may then be calculated as CI - DD.
%
% If the observed process covariance matrix G and/or per-element state prediction
% error covariance matrices P are empty, they are calculated. If many projections L
% are to be processed, precomputing G and P speeds up computation considerably. See
% iss2ce_precomp.m.
%
% L        orthonormal subspace basis
% A,C,K,V  ISS parameters
% G        ISS covariance matrix (if empty, calculate)
% P        per-element state prediction error covariance matrices (if empty, calculate)
%
% NOTE: does NOT assume normalisation of residuals covariance matrix V.

noG = nargin < 6 || isempty(G);
noP = nargin < 7 || isempty(P);

[n,m] = size(C);

VCHOL = chol(V,'lower');
KV  = K*VCHOL;
KVK = KV*KV';
LV  = L'*VCHOL;
LVL = LV*LV';

% First term

if noP
	P = zeros(m,m,n);
	for i = 1:n
		[~,~,~,~,P(:,:,i)] = mdare(A,C(i,:),KVK,V(i,i),K*V(:,i));
	end
end
i1 = zeros(n,1);
for i = 1:n
	i1(i) = logdet(L'*(C*P(:,:,i)*C'+V)*L);
end
I1 = sum(i1);

% Second term

[~,VR] = mdare(A,L'*C,KVK,LVL,KV*LV');
I2 = logdet(VR);

% Third term

if noG
	M = dlyap(A,KVK); % state covariance
	G  = C*M*C' + V;  % observable covariance
end
I3 = (n-1)*logdet(L'*G*L);

% Fourth term

I4 = logdet(LVL);

% Co-information

CI = I1 - I4 - I3;

% Dynamical Dependence

DD = I2 - I4;

% Causal Emergence

CE = I1 - I2 - I3; % = CI - DD
