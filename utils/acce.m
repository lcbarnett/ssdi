function CE = acce(L,G,A,C,K)

% Calculate causal emergence of projection L from autocovariance
% sequence G.

[n,n1,p1] = size(G);
assert(n1 == n);
p = p1-1;

[n1,m] = size(L);
assert(n1 == n);

G0yy = L'*G(:,:,1)*L;

Gyx = zeros(p,m);
Hyyx = zeros(n,1);
for i = 1:n

	for k = 1:p
		Gyx(k,:) = L'*G(:,i,k+1);
	end
	Gxx = G(i,i,1:p);
	Ayx = tsolve(Gxx(:)',Gyx)';
	Vyyx = G0yy-Ayx*Gyx;
	Hyyx(i) = logd(Vyyx);

end

[~,Vyy,rep] = ssdare(A,L'*C,K*K',[],K*L);
if rep < 0  || rep > 1e-08 % DARE failed
	Hyy = NaN;
else
	Hyy = logd(Vyy);
end

CE = Hyy-sum(Hyyx);

function LD = logd(V)

[L,flag] = chol(V);
if flag == 0 % symmetric, positive-definite
	LD = 2*sum(log(diag(L)));
else
	LD = NaN;
end
