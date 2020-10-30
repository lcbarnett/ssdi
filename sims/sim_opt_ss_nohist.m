%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('resdir', 'var'), resdir = tempdir; end % results directory
if ~exist('rid',    'var'), rid    = '';      end % run ID tag

if exist('G','var')
	varmod = true;
	n = size(G,1); % microscopic state dimension
	if ~exist('r','var'), r = n; end % var model order
	if ~exist('w','var'), w = 1; end % AR coefficients decay parameter
else
	varmod = false;
	if ~exist('n','var'), n = 7;   end % microscopic state dimension
	if ~exist('r','var'), r = 3*n; end % hidden state dimension
	G = ones(n); % fully-connected
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

scriptname = mfilename;

if es1p1, algo = '1+1 ES'; else, algo = 'SD'; end

% Generate random model

mrstate = rng_seed(mseed);
V = eye(n); % residuals covariance is decorrelated and normalised to unity
if varmod
	ARA = var_rand(G,r,rho,w);
	[A,C,K] = var_to_ss(ARA);
	gc = var_to_pwcgc(ARA,V);
else
	[A,C,K] = iss_rand(n,r,rho);
	gc = ss_to_pwcgc(A,C,K,V);
end
rng_restore(mrstate);

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

orstate = rng_seed(oseed);
for k = 1:nruns
	fprintf('run %2d of %2d ... ',k,nruns);
	if es1p1
		[dopt(k),converged,sig,iopt(k),Lopt(:,:,k)] = opt_ss_es1p1(A,C,K,L0(:,:,k),niters,sig0,ifac,nfac,estol);
		fprintf('dopt = %.4e : sig = %.4e : ',dopt(k),sig);
		if converged, fprintf('converged in %d iterations\n',iopt(k)); else, fprintf('unconverged\n'); end
	else
		[dopt(k),Lopt(:,:,k)] = opt_ss_sd(A,C,K,L0(:,:,k),niters,sig0);
		fprintf('dopt = %.4e\n',dopt(k));
	end

end
rng_restore(orstate);

% Sort (local) optima by dynamical dependence

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
fprintf('\noptimal dynamical dependence =\n'); disp(dopt');

clear k

wsfile = fullfile(resdir,[scriptname rid '.mat']);
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');
