%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('resdir', 'var'), resdir = tempdir; end % results directory
if ~exist('rid',    'var'), rid    = '';      end % run ID tag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('r',       'var'), r       = 5;       end % AR lags
if ~exist('rho',     'var'), rho     = 0.9;     end % spectral norm (< 1)
if ~exist('mseed',   'var'), mseed   = 0;       end % model random seed (0 to use current rng state)
if ~exist('sig0',    'var'), sig0    = 0.1;     end % (initial) step size
if ~exist('es1p1',   'var'), es1p1   = true;    end % use 1+1 evolutionary strategy? (else stochastic gradient descent)
if ~exist('esrule',  'var'), esrule  = 1/5;     end % evolutionary strategy adaptation rule
if ~exist('estol',   'var'), estol   = 1e-8;    end % evolutionary strategy convergence tolerance
if ~exist('niters',  'var'), niters  = 1000;    end % iterations
if ~exist('nruns',   'var'), nruns   = 10;      end % runs (restarts)
if ~exist('iseed',   'var'), iseed   = 0;       end % initialisation random seed (0 to use current rng state)
if ~exist('oseed',   'var'), oseed   = 0;       end % optimisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',  'var'), gvprog  = 'fdp';   end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',  'var'), gvdisp  = true;    end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',  'var'), gpterm  = 'x-pdf'; end % Gnuplot terminal
if ~exist('gpscale', 'var'), gpscale = 1.2;     end % Gnuplot scale factor(s)
if ~exist('gpfsize', 'var'), gpfsize = 14;      end % Gnuplot font size
if ~exist('gpplot',  'var'), gpplot  = 2;       end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;

if es1p1, algo = '1+1 ES'; else, algo = 'SD'; end

% Generate random model

n = 7;
nx = 2;
ny = 3;
nz = n-nx-ny;
mrstate = rng_seed(mseed);
x = 1:nx;
y = nx+1:nx+ny;
z = nx+ny+1:n;
ARA = randn(n,n,r);
ARA(x,[y z],:) = zeros(nx,ny+nz,r);
ARA
T = orth(randn(n));
for p = 1:r, ARA(:,:,p) = T*ARA(:,:,p)*T'; end
ARA(y,[x z],:) = zeros(ny,nx+nz,r);
ARA
for p = 1:r, ARA(:,:,p) = T'*ARA(:,:,p)*T; end
ARA
return
ARA = specnorm(ARA,rho);
[A,C,K] = var_to_ss(ARA);
V = eye(n); % residuals covariance is decorrelated and normalised to unity
gc = var_to_pwcgc(ARA,V);
rng_restore(mrstate);


% Calculate causal graph, and display

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);

%return

% Set 1+1 evolutionary strategy parameters

if es1p1
	[ifac,nfac] = es1p1_facs(esrule,m*(n-m));
end

% Initial linear projections

irstate = rng_seed(iseed);
L0 = zeros(n,m,nruns);
for k = 1:nruns
	L0(:,:,k) = orthonormalise(randn(n,m));
end
rng_restore(irstate);

iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = zeros(n,m,nruns);

dd = nan(niters,nruns);

orstate = rng_seed(oseed);
for k = 1:nruns
	fprintf('run %2d of %2d ... ',k,nruns);
	if es1p1
		[dopt(k),dd(:,k),converged,sig,iopt(k),Lopt(:,:,k)] = opt_ss_es1p1_x(A,C,K,L0(:,:,k),niters,sig0,ifac,nfac,estol);
		fprintf('dopt = %.4e : sig = %.4e : ',dopt(k),sig);
		if converged, fprintf('converged in %d iterations\n',iopt(k)); else, fprintf('unconverged\n'); end
	else
		[dopt(k),dd(:,k),Lopt(:,:,k)] = opt_ss_sd_x(A,C,K,L0(:,:,k),niters,sig0);
		fprintf('dopt = %.4e\n',dopt(k));
	end

end
rng_restore(orstate);

% Sort (local) optima by dynamical dependence

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
fprintf('\noptimal dynamical dependence =\n'); disp(dopt');

nweight = zeros(n,nruns);
for k = 1:nruns
	nweight(:,k) = 1-gmetricsx(Lopt(:,:,k));
end

Loptd = gmetrics(Lopt);

clear k

wsfile = fullfile(resdir,[scriptname rid '.mat']);
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');

% Plot optimisation histories

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_opthist' rid]);
gp_opthist(dd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
