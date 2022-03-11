function [uidx,usiz,nruns] = Lcluster(dist,tol,dd,gpterm,gpscale,gpfsize,gpplot,logsy,resdir,rid);

% Hyperplanes should be sorted (ascending) by
% dynamical dependence prior to calling.

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
		uopt(k) = dd(i);
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
	fprintf('run %4d      size = %3d      dd = %7.4f\n',uidx(k),usiz(k),uopt(k));
end
fprintf('------------------------------------------\n\n');

if nargin > 3 && ~isempty(dd) && ~isempty(gpterm)
	if nargin <  5, gpscale = []; end
	if nargin <  6, gpfsize = []; end
	if nargin <  7, gpplot  = []; end
	if nargin <  8 || isempty(logsy),  logsy  = true;    end
	if nargin <  9 || isempty(resdir), resdir = tempdir; end
	if nargin < 10 || isempty(rid),    rid    = '';      end
	gpdata = [(1:r)' diff([0;dd'])];
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
