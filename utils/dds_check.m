function derr = dds_check(A,C,K,H,m,nsics,verb)

% Check accuracy of spectral DD against state-space DD
% at given scale m across a set of random projections.

if nargin < 7 || isempty(verb), verb = 1; end

fres = size(H,3)-1;
if verb > 0
	fprintf('Spectral DD accuracy check (m = %d, frequency resolution = %d) ... ',m,fres);
end

n = size(C,1);
L = rand_orthonormal(n,m,nsics); % (orthonormalised) random linear projections
D1 = zeros(nsics,1);
for k = 1:nsics
	D1(k) = iss2dd(L(:,:,k),A,C,K);
end
D2 = zeros(nsics,1);
for k = 1:nsics
	D2(k) = trfun2dd(L(:,:,k),H);
end
derr = maxabs(D1-D2);

if verb > 0
	fprintf('integration error = %e\n',derr);
	if derr > 1e-12, fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n'); end
	fprintf('\n');
end
