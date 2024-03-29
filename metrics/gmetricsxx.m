function d = gmetricsxx(L,maxangle)

if nargin < 2 maxangle = []; end % gmetric default

[n,m] = size(L);

c = nchoosek(1:n,m); % combinations
nc = size(c,1);      % number of combinations
d = zeros(nc,1);
for k = 1:nc
	v = zeros(n,m);
	v(c(k,:),:) = eye(m);
	d(k) = gmetric(L,v,maxangle);
end
