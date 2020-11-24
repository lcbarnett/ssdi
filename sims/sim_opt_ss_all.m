%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation - all scales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if ~exist('gvprog',  'var'), gvprog  = 'dot';   end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',  'var'), gvdisp  = true;    end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',  'var'), gpterm  = 'x-pdf'; end % Gnuplot terminal
if ~exist('gpscale', 'var'), gpscale = 1.2;     end % Gnuplot scale factor(s)
if ~exist('gpfsize', 'var'), gpfsize = 14;      end % Gnuplot font size
if ~exist('gpplot',  'var'), gpplot  = 2;       end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;

if es1p1, algo = '1+1 ES'; else, algo = 'SD'; end

n1 = n-1;

% Generate random model

mrstate = rng_seed(mseed);
V = eye(n); % residuals covariance is decorrelated and normalised to unity
if varmod
	ARA = var_rand(G,r,rho,w);
%{
T = orth(randn(n));
for u=1:n
	ARA(:,:,u) = T*ARA(:,:,u)*T';
	ARA(3:7,1:2,u) = zeros(5,2);
	ARA = specnorm(ARA,rho);
end
%}
	[A,C,K] = var_to_ss(ARA);
	gc = var_to_pwcgc(ARA,V);
else
	[A,C,K] = iss_rand(n,r,rho);
	gc = ss_to_pwcgc(A,C,K,V);
end
rng_restore(mrstate);

% Calculate causal graph, and display

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);

% Initial linear projections

irstate = rng_seed(iseed);
P0 = cell(n1,1);
for m = 1:n1
	P0{m} = randn(n,m,nruns);
end
rng_restore(irstate);

iopt = zeros(nruns,n1);
dopt = zeros(nruns,n1);
Lopt = cell(n1,1);
for m = 1:n1
	Lopt{m} = zeros(n,m,nruns);
end

if es1p1 % Set 1+1 evolutionary strategy parameters
	ifac = zeros(n1,1);
	nfac = zeros(n1,1);
	for m = 1:n1
		[ifac(m),nfac(m)] = es1p1_facs(esrule,m*(n-m));
	end
end

orstate = rng_seed(oseed);
% parfor m = 1:n1 % for some reason, MATLAB doesn't complain, but then doesn't parallelise this loop!
for m = 1:n1
	fprintf('m = %2d of %2d\n',m,n1);
	for k = 1:nruns
		fprintf('\trun %2d of %2d ... ',k,nruns);
		if es1p1
			[dopt(k,m),converged,sig,iopt(k,m),Lopt{m}(:,:,k)] = opt_ss_es1p1(A,C,K,P0{m}(:,:,k),niters,sig0,ifac(m),nfac(m),estol);
			fprintf('dopt = %.4e : sig = %.4e : ',dopt(k,m),sig);
			if converged, fprintf('converged in %d iterations\n',iopt(k,m)); else, fprintf('unconverged\n'); end
		else
			[dopt(k,m),Lopt{m}(:,:,k)] = opt_ss_sd(A,C,K,P0{m}(:,:,k),niters,sig0);
			fprintf('dopt = %.4e\n',dopt(k,m));
		end
	end
end
rng_restore(orstate);

% Sort (local) optima by dynamical dependence

for m = 1:n1
	[dopt(:,m),sidx] = sort(dopt(:,m));
	iopt(:,m) = iopt(sidx,m);
	Lopt{m} = Lopt{m}(:,:,sidx);
end

fprintf('\noptimal dynamical dependence =\n\n');
disp(num2str(1:n1,'    %8d'));
disp(num2str(dopt,'    %8.6f'));

eweight = gc/nanmax(gc(:));

nweight = zeros(n,n1,nruns);
for m = 1:n1
	for k = 1:nruns
		nweight(:,m,k) = 1-gmetricsx(Lopt{m}(:,:,k));
	end
end

Loptd = zeros(nruns,nruns,n1);
for m = 1:n1
	Loptd(:,:,m) = gmetrics(Lopt{m});
end

clear m k

wsfile = fullfile(resdir,[scriptname rid '.mat']);
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');

% Plot (local) optimum dynamical dependencies at all scales

gptitle = sprintf('Local optima (%s) : n = %d, r = %d',algo,n,r);
gpstem = fullfile(resdir,[scriptname '_localopt' rid]);
gp_localopt(dopt,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% m = 4; k = 3; wgraph2dot(nweight(:,m,k),eweight,fullfile(resdir,sprintf('graph%s_scale%02d_run%03d',rid,m,k)),[],gvprog);
