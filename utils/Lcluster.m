function [uidx,usiz] = Lcluster(dist,tol,dd,logsy,resdir,rid,gpterm,gpscale,gpfsize,gpplot);

% Hyperplanes should be sorted (ascending) by
% dynamical dependence prior to calling.

r = size(dist,1);
a = true(1,r); % still available
uidx = zeros(1,r);
usiz = zeros(1,r);
k = 0;
for i = 1:r
	if a(i) % new cluster
		k = k+1;
		uidx(k) = i;
		usiz(k) = 1;
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

if nargin > 2 && ~isempty(dd)
	if nargin <  4 || isempty(logsy),  logsy  = true;    end
	if nargin <  5 || isempty(resdir), resdir = tempdir; end
	if nargin <  6 || isempty(rid),    rid    = '';      end
	if nargin <  7, gpterm  = []; end
	if nargin <  8, gpscale = []; end
	if nargin <  9, gpfsize = []; end
	if nargin < 10, gpplot  = []; end
	gpdata = [(1:r)' diff([0;dd'])];
	gpname  = ['lcluster' rid];
	gpstem  = fullfile(resdir,gpname);
	gp_write(gpstem,gpdata);
	gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
	fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
	fprintf(gp,'set title "%s"\n','DD difference');
	fprintf(gp,'unset key\n');
	fprintf(gp,'set grid\n');
	fprintf(gp,'set xlabel "run number"\n');
	if logsy
		fprintf(gp,'set logs y\n');
	end
	fprintf(gp,'set style fill solid 0.2\n');
	fprintf(gp,'\nplot datfile using 1:2 with boxes\\\n');
	gp_close(gp,gpstem,gpterm,gpplot);
end
