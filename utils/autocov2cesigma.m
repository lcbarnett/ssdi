function Vi = autocov2cesigma(Gamma,i)

% Supply autocovariance sequence G and variable index i

Gi  = squeeze(Gamma(:,i,2:end));
Gii = toeplitz(squeeze(Gamma(i,i,1:end-1)));
Fi  = Gi/chol(Gii,'lower');
G0c = chol(Gamma(:,:,1),'lower');
Vi  = G0c*G0c'-Fi*Fi';
