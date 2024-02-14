function [CESRC,CRC] = ac2ces(G)

% Calculate the CE Sigma_i matrices
%
% G        autocovariance sequence
% CRC      right (upper) Cholesky factor of covariance matrix
% CESRC    right (upper) Cholesky factors of the CE Sigma_i matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: all calculations in UNTRANSFORMED coordinates!!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CRC = chol(G(:,:,1));
n = size(G,1);
CESRC = zeros(n,n,n);
for i = 1:n
	fprintf('.');
	Gi  = squeeze(G(:,i,2:end));
	Gii = toeplitz(squeeze(G(i,i,1:end-1)));
	Fi  = Gi/chol(Gii);
	CESRC(:,:,i) = chol(CRC'*CRC-Fi*Fi');
end
