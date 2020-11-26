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
if ~exist('piters',  'var'), piters  = 10000;   end % pre-optimisation iterations
if ~exist('nruns',   'var'), nruns   = 10;      end % runs (restarts)
if ~exist('nnorm',   'var'), nnorm   = 1000;    end % normalisation projections
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

% Calculate causal graph, and display

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);

%return

% Set 1+1 evolutionary strategy parameters

if es1p1
	[ifac,nfac] = es1p1_facs(esrule,m*(n-m));
end

irstate = rng_seed(iseed);

% Normalisation factors

if nnorm > 0
	L = randn(n,m,nnorm);
	D = zeros(nnorm,1);
	if piters > 0
		CAK = CAK_seq(A,C,K);
		for k = 1:nnorm
			[Lk,Mk] = orthonormalise(L(:,:,k));
			D(k) = ssdd1(Lk,Mk,CAK);
		end
		ddpmean = mean(D);
		fprintf('"fake" dd mean = %g : ',ddpmean);
	end
	for k = 1:nnorm
		Lk = orthonormalise(L(:,:,k));
		D(k) = ssdd(Lk,A,C,K);
	end
	ddnmean = mean(D);
	fprintf('dd mean = %g\n\n',ddnmean);
end

% Initial linear projections

P0 = randn(n,m,nruns);

rng_restore(irstate);

iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = zeros(n,m,nruns);

titers = piters+niters;
dd = nan(titers,nruns);

orstate = rng_seed(oseed);
for k = 1:nruns
	fprintf('run %2d of %2d ... ',k,nruns);
	if es1p1
		if piters > 0
			[doptk,ddpk,converged,sigk,ioptk,P0(:,:,k)] = opt_ss1_es1p1_x(A,C,K,P0(:,:,k),piters,sig0,ifac,nfac,estol);
			fprintf('dopt = %.4e : sig = %.4e : ',doptk,sigk);
			if converged, fprintf('preopt converged in %d iterations : ',ioptk); else, fprintf('preopt unconverged : '); end
			if nnorm > 0, ddpk = ddpk/ddpmean; end
		end
		[dopt(k),ddnk,converged,sig,iopt(k),Lopt(:,:,k)] = opt_ss_es1p1_x(A,C,K,P0(:,:,k),niters,sig0,ifac,nfac,estol);
		fprintf('dopt = %.4e : sig = %.4e : ',dopt(k),sig);
		if converged, fprintf('converged in %d iterations\n',iopt(k)); else, fprintf('unconverged\n'); end
		if nnorm > 0, ddnk = ddnk/ddnmean; end
	else
		if piters > 0
			[doptk,ddpk,converged,sigk,ioptk,P0(:,:,k)] = opt_ss1_sd_x(A,C,K,P0(:,:,k),piters,sig0);
			fprintf('dopt = %.4e : sig = %.4e : ',doptk,sigk);
			if converged, fprintf('preopt converged in %d iterations : ',ioptk); else, fprintf('preopt unconverged : '); end
			if nnorm > 0, ddpk = ddpk/ddpmean; end
		end
		[dopt(k),ddnk,Lopt(:,:,k)] = opt_ss_sd_x(A,C,K,P0(:,:,k),niters,sig0);
		fprintf('dopt = %.4e\n',dopt(k));
		if nnorm > 0, ddnk = ddnk/ddnmean; end
	end
	ddk = [ddpk(~isnan(ddpk));ddnk(~isnan(ddnk))];
	dd(1:length(ddk),k) = ddk;
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
