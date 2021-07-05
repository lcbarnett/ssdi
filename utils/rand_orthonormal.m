function [L,M] = rand_orthonormal(n,m,r)

% Returns r random n x m orthonormal matrices, and
% optionally their n x (n-m) orthogonal complements.

L = randn(n,m,r);
if nargout > 1
	M = zeros(n,m,r);
	for k = 1:r
		[L(:,:,k),M(:,:,k)] = orthonormalise(L(:,:,k));
	end
else
	for k = 1:r
		L(:,:,k) = orthonormalise(L(:,:,k));
	end
end
