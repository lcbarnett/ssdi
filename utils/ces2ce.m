function  CE = ces2ce(L,CES,G0c,DD)

% Calculate causal emergence of projection L
%
% To calculate CE sequence CES and covariance matrix
% lower Cholesky factor G0c, see ac2ces.
%
% NOTE 1: assumes uncorrelated residuals
% NOTE 2: projection L MUST be orthonormal!!!

n = size(G0c,1);
CEH = zeros(n,1);
for i = 1:n
	CEH(i) = logdet(L'*CES(:,:,i)*L);
end
LG0c = L'*G0c;
CE = -(n-1)*logdet(LG0c*LG0c') - DD + sum(CEH);
