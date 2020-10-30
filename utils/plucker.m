function p = plucker(L,nrmlz)

% Plucker embedding

[n,m] = size(L);

c = nchoosek(1:n,m); % combinations
nc = size(c,1);      % number of combinations
p = zeros(nc,1);
for k = 1:nc
	p(k) = det(L(c(k,:),:));
end

if nargin < 2 || nrmlz
	p = p/sqrt(sum(p.^2));
end
