%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default parameters (override on command line)

if ~exist('n',       'var'), n       = 7;       end % microscopic state dimension
if ~exist('r',       'var'), r       = 3*n;     end % hidden state dimension
if ~exist('rhoa',    'var'), rhoa    = 0.9;     end % state-space AR spectral norm (< 1)
if ~exist('m',       'var'), m       = 4;       end % macroscopic state dimension
if ~exist('sig0',    'var'), sig0    = 0.1;     end % (initial) step size
if ~exist('es1p1',   'var'), es1p1   = true;    end % use 1+1 evolutionary strategy? (else stochastic gradient descent)
if ~exist('esrule',  'var'), esrule  = 1/5;     end % evolutionary strategy adaptation esrule
if ~exist('esdim',   'var'), esdim   = m*(n-m); end % evolutionary strategy effective search-space dimension
if ~exist('estol',   'var'), estol   = 1e-8;    end % evolutionary strategy convergence tolerance
if ~exist('iters',   'var'), iters   = 1000;    end % iterations
if ~exist('runs',    'var'), runs    = 10;      end % runs (restarts)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gpterm',  'var'), gpterm  = 'x-pdf'; end % Gnuplot terminal
if ~exist('gpscale', 'var'), gpscale = 1.2;     end % Gnuplot scale factor(s)
if ~exist('gpfsize', 'var'), gpfsize = 14;      end % Gnuplot font size

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if es1p1
	% Set 1+1 evolutionary strategy parameters
	if esrule < eps  % special case (ignore gain): esrule = 1/5, gain = (5/4)*log(3/2) = 0.5068...
		ifac = 3/2;
		nfac = realpow(3/2,-1/4);
	else
		gain = 1/sqrt(esdim);
		ifac = exp((1-esrule)*gain);
		nfac = exp(-esrule*gain);
	end
end

% Generate a random state-space model in innovations form

[A,C,K,rhob] = iss_rand(n,r,rhoa);

% Initial linear projections

L0 = zeros(n,m,runs);
for k = 1:runs
	L0(:,:,k) = orthonormalise(randn(n,m));
end

iopt = zeros(1,runs);
dopt = zeros(1,runs);
Lopt = zeros(n,m,runs);

dd = nan(iters,runs);

for k = 1:runs

	fprintf('run %2d of %2d ... ',k,runs);

	% Initialise

	L = L0(:,:,k);

	% Calculate dynamical dependence of initial projection

	d = ssdd(L,A,C,K);

	dd(1,k) = d;

	sig = sig0;
	converged = false;

	% Optimise

	for i = 2:iters

		% "Mutate" projection and orthonormalise
		Ltry = orthonormalise(L + sig*randn(n,m));

		% Calculate dynamical dependence of mutated projection

		dtry = ssdd(Ltry,A,C,K);

		% If dynamical dependence smaller, accept mutant

		if es1p1
			if dtry < d
				L = Ltry;
				d = dtry;
				sig = ifac*sig;
			else
				sig = nfac*sig;
			end
		else
			if dtry < d
				L = Ltry;
				d = dtry;
			end
		end
		dd(i,k) = d;

		if sig < estol
			converged = true;
			break
		end

	end

	dopt(k)     = d;
	iopt(k)     = i;
	Lopt(:,:,k) = L;

	if es1p1
		if converged, fprintf('converged in %d iterations\n',iopt(k)); else, fprintf('unconverged (sig = %g)\n',sig); end
	else
		fprintf('done\n');
	end
end

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
fprintf('\noptimal dynamical dependence =\n'); disp(dopt');

ii = (1:iters)';
ddmax = max(dd(:));
if es1p1, algo = '1+1 ES'; else, algo = 'SD'; end

gpstem = fullfile(tempdir,mfilename);
gp_write(gpstem,[ii dd]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpstem);
fprintf(gp,'set title "%s : n = %d, r = %d, m = %d"\n',algo,n,r,m);
fprintf(gp,'set key top left Left rev\n');
fprintf(gp,'set xr [1:%g]\n',iters);
fprintf(gp,'set yr [0:%g]\n',1.05*ddmax);
fprintf(gp,'set xlabel "iterations"\n');
fprintf(gp,'set ylabel "dynamical dependence"\n');
fprintf(gp,'set logs x\n');
fprintf(gp,'set grid\n');
fprintf(gp,'plot \\\n');
for r = 1:runs
	fprintf(gp,'datfile using 1:%d w lines not, \\\n',r+1);
end
fprintf(gp,'NaN not\n');
gp_close(gp,gpstem,gpterm);
