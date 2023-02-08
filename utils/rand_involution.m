function Q = rand_involution(n,m,r)

% Returns r random n x n involution matrices with trace = 2*m-n
%
% Ref: Z. Lai, L.-H. Lim and  K.Ye, "Simpler Grassmannian optimization"
% arXiv:2009.13502 [math.OC] - https://arxiv.org/abs/2009.13502

if nargin < 3 || isempty(r), r = 1; end

J = eye(n);
J(m+1:n,m+1:n) = -eye(n-m);
Q = zeros(n,n,r);
for k = 1:r
	V = orthonormalise(randn(n));
	Q(:,:,k) = V*J*V';
end
