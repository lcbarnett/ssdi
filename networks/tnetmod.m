function G = tnetmod(smod,cons)

nmods = length(smod);
nnodes = sum(smod);
G = eye(nnodes);

mods = cell(nmods,1);
k = 0;
for m = 1:nmods
	mods{m} = k+1:k+smod(m);
	k = k+smod(m);
end

% Fully intra-connect modules
for m = 1:nmods
	G(mods{m},mods{m}) = ones(smod(m));
end

[ncons,w] = size(cons);
if w == 2 % fully inter-connect specified modules
	for c = 1:ncons
		G(mods{cons(c,1)},mods{cons(c,2)}) = ones(smod(cons(c,1)),smod(cons(c,2)));
	end
else      % inter-module connections
	for c = 1:ncons
		G(mods{cons(c,1)}(cons(c,2)),mods{cons(c,3)}(cons(c,4))) = 1;
	end
end
