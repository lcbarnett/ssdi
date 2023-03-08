function [CES,G0c] = ac2ces(G)

% Supply autocovariance sequence G.

G0c = chol(G(:,:,1),'lower');
n = size(G,1);
CES = zeros(n,n,n);
for i = 1:n
	fprintf('.');
	Gi  = squeeze(G(:,i,2:end));
	Gii = toeplitz(squeeze(G(i,i,1:end-1)));
	Fi  = Gi/chol(Gii,'lower');
	CES(:,:,i) = G0c*G0c'-Fi*Fi';
end
