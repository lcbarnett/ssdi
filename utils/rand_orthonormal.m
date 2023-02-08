function [L,M] = rand_orthonormal(n,m,r)

% Returns r random n x m orthonormal matrices, and
% optionally their n x (n-m) orthogonal complements.

if nargin < 3 || isempty(r), r = 1; end

L = randn(n,m,r);
if nargout > 1
	M = zeros(n,n-m,r);
	for k = 1:r
		[L(:,:,k),M(:,:,k)] = orthonormalise(L(:,:,k));
	end
else
	for k = 1:r
		L(:,:,k) = orthonormalise(L(:,:,k));
	end
end
