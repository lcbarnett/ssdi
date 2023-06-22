function d = gmetrics(L,maxangle)

if nargin < 2 maxangle = []; end % gmetric default

N = size(L,3);
d = zeros(N);
for r1 = 1:N;
	for r2 = r1+1:N
		d(r1,r2) = gmetric(L(:,:,r1),L(:,:,r2),maxangle);
	end
end
d = symmetrise(d);
