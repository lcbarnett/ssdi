function  [CE,DD,COI] = ces2ce(L,H,VRC,CESRC,CRC)

% Calculate causal emergence of projection L
%
% L        orthonormal subspace basis
% H        transfer function
% VRC      right (upper) Cholesky factor of residuals covariance matrix
% CESRC    right (upper) Cholesky factors of the CE Sigma_i matrices, as returned by ac2ces
% CRC      right (upper) Cholesky factor of covariance matrix, as returned by ac2ces
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: all calculations in UNTRANSFORMED coordinates!!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate reduced residuals generalised covariance

h = size(H,3);
v = zeros(h,1);
for k = 1:h % over [0,pi]
	RHLk = VRC*H(:,:,k)'*L;
	v(k) = sum(log(diag(chol(RHLk'*RHLk)))); % (log-determinant)/2
end
VRED = sum(v(1:end-1)+v(2:end))/(h-1); % integrate frequency-domain DD (trapezoidal rule) to get time-domain DD

% Calculate causal emergence

n = size(CRC,1);
CEH = zeros(n,1);
for i = 1:n
	CESRCLi = CESRC(:,:,i)*L;
	CEH(i) = logdet(CESRCLi'*CESRCLi);
end
CRCL = CRC*L;
CE = -(n-1)*logdet(CRCL'*CRCL) + sum(CEH) - VRED;

% Optionally return dynamical dependence

if nargout > 1
	RHL = VRC*L;
	DD = VRED - 2*sum(log(diag(chol(RHL'*RHL))));
	if nargout > 2
		COI = CE+DD;
	end
end
