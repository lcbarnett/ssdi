function d = gmetricsx(L)

n = size(L,1);
d = zeros(n,1);
for i = 1:n
	v = zeros(n,1);
	v(i) = 1;
	d(i) = gmetric(L,v);
end
