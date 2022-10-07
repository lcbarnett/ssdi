function CES = ac2ces(Gamma)

% Supply autocovariance sequence Gamma.

G0c = chol(Gamma(:,:,1),'lower');
n = size(Gamma,1);
CES = zeros(n,n,n);
for i = 1:n
	fprintf('.');
	Gi  = squeeze(Gamma(:,i,2:end));
	Gii = toeplitz(squeeze(Gamma(i,i,1:end-1)));
	Fi  = Gi/chol(Gii,'lower');
	CES(:,:,i) = G0c*G0c'-Fi*Fi';
end

% To calculate CE using DARE for DD
%
% CEH = zeros(n,1);
% for i = 1:n
% 	CEH(i) = logdet(L'*CES(:,:,i)*L);
% end
% DD = iss2dd(L,A,C,K);
% LG0 = L'*chol(Gamma(:,:,1),'lower');
% CE = -(n-1)*logdet(LG0*LG0') - DD + sum(CEH);
%
% To calculate CE using spectral method for DD, use
%
% DD = trfun2dd(L,H);
