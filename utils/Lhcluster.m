function [clust,cddep,csize,T,Z] = Lhcluster(L,ddep,coff,method,mangle,fignum);

% Hierarchical clustering of hyperplanes by distance (angle metric)

if nargin < 4 || isempty(method), method = 'average'; end
if nargin < 5 || isempty(mangle), mangle = true;      end
if nargin < 6 || isempty(fignum), fignum = 1;         end

[n,m,N] = size(L);
LL = zeros(N,n*m);
for k = 1:N
	Lk = L(:,:,k);
	LL(k,:) = Lk(:);
end

% Create

Z = linkage(LL,method,@distfun);

figure(fignum); clf
dendrogram(Z);

T = cluster(Z,'Cutoff',coff);

nclust = max(T);
clust = cell(nclust,1);
cddep = cell(nclust,1);
csize = zeros(nclust,1);
for c = 1:nclust
	clust{c} = find(T == c)';
	csize(c) = numel(clust{c});
	cdd = ddep(clust{c});
	[~,didx] = sort(cdd,'ascend'); % sort by ascending DD within cluster
	clust{c} = clust{c}(didx);
	cddep{c} = cdd(didx);
end
[~,sidx] = sort(csize,'descend'); % sort clusters by descending size
clust = clust(sidx);
cddep = cddep(sidx);
csize = csize(sidx);



%{
r = size(dist,1);
a = true(1,r); % still available
uidx = zeros(1,r);
usiz = zeros(1,r);
uopt = zeros(1,r);
k = 0;
for i = 1:r
	if a(i) % new cluster
		k = k+1;
		uidx(k) = i;
		usiz(k) = 1;
		uopt(k) = ddep(i);
		a(i) = false;
		for j = 1:r
			if a(j)
				if dist(i,j) < tol % add j to cluster
					a(j) = false;
					usiz(k) = usiz(k)+1;
				end
			end
		end
	end
end
uidx = uidx(1:k);
usiz = usiz(1:k);
uopt = uopt(1:k);

nruns = length(uidx);
fprintf('\n------------------------------------------\n');
fprintf('Clusters = %d\n',nruns);
fprintf('------------------------------------------\n');
for k = 1:nruns
	fprintf('run %4d      size = %3d      ddep = %7.4f\n',uidx(k),usiz(k),uopt(k));
end
fprintf('------------------------------------------\n\n');

if nargin > 3 && ~isempty(ddep) && ~isempty(gpterm)
	if nargin <  5, gpscale = []; end
	if nargin <  6, gpfsize = []; end
	if nargin <  7, gpplot  = []; end
	if nargin <  8 || isempty(logsy),  logsy  = true;    end
	if nargin <  9 || isempty(resdir), resdir = tempdir; end
	if nargin < 10 || isempty(rid),    rid    = '';      end
	gpdata = [(1:r)' diff([0;ddep'])];
	gpname  = ['lcluster' rid];
	gpstem  = fullfile(resdir,gpname);
	gp_write(gpstem,gpdata);
	gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
	fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
	fprintf(gp,'set title "%s"\n','DD difference');
	fprintf(gp,'unset key\n');
	fprintf(gp,'set grid\n');
	fprintf(gp,'set xr [0:%g]\n',r+1);
	fprintf(gp,'set xlabel "run number"\n');
	if logsy
		fprintf(gp,'set logs y\n');
	end
	fprintf(gp,'set style fill solid 0.2\n');
	fprintf(gp,'\nplot datfile using 1:2 with boxes\\\n');
	gp_close(gp,gpstem,gpterm,gpplot);
	fprintf('\n');
end
%}

	function d = distfun(L1,L2) % Metric on the Grassmanian manifold, normalised to lie in [0,1]
		L1 = reshape(L1,n,m);
		k = size(L2,1);
		d = zeros(k,1);
		if mangle
			for r = 1:k
				L2r = reshape(L2(r,:),n,m);
				theta = subspacea(L1,L2r);
				d(r) = max(theta)/(pi/2);
			end
		else
			for r = 1:k
				L2r = reshape(L2(r,:),n,m);
				theta = subspacea(L1,L2r);
				d(r) = sqrt(mean(theta.^2))/(pi/2);
			end
		end
	end

end
