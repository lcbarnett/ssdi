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

if ~exist('rho',      'var'), rho      = 0.9;       end % spectral norm (< 1)
if ~exist('mseed',    'var'), mseed    = 0;         end % model random seed (0 to use current rng state)
if ~exist('sig0',     'var'), sig0     = 0.1;       end % (initial) step size
if ~exist('esrule',   'var'), esrule   = 1/5;       end % evolutionary strategy adaptation rule
if ~exist('estol',    'var'), estol    = 1e-8;      end % evolutionary strategy convergence tolerance
if ~exist('piters',   'var'), piters   = 10000;     end % pre-optimisation 1 iterations
if ~exist('siters',   'var'), siters   = 1000;      end % pre-optimisation 2 iterations
if ~exist('niters',   'var'), niters   = 100;       end % iterations
if ~exist('nruns',    'var'), nruns    = 10;        end % runs (restarts)
if ~exist('sigreset', 'var'), sigreset = false;     end % reset step size between optimisations?
if ~exist('nnorm',    'var'), nnorm    = 1000;      end % normalisation projections
if ~exist('ddscheck', 'var'), ddscheck = 100;       end % iterations for spectral integration check
if ~exist('hist',     'var'), hist     = true;      end % calculate optimisation history?
if ~exist('iseed',    'var'), iseed    = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('oseed',    'var'), oseed    = 0;         end % optimisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',   'var'), gvprog   = 'dot';     end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',   'var'), gvdisp   = true;      end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpscale',  'var'), gpscale  = [Inf 0.8]; end % Gnuplot scale factor(s)
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;

% Generate random VAR or ISS model

mrstate = rng_seed(mseed);
V = eye(n); % residuals covariance is decorrelated and normalised to unity
if varmod
	ARA = var_rand(G,r,rho,w);
	[A,C,K] = var_to_ss(ARA);
	CAK = ARA;                  % CAK sequence for pre-optimisation
	info = var_info(ARA,V);     % VAR information
	fres = info.fres;           % frequency resolution
	H = var2trfun(ARA,fres);    % transfer function
	gc = var_to_pwcgc(ARA,V);   % causal graph
else
	[A,C,K] = iss_rand(n,r,rho);
	CAK = iss2cak(A,C,K);       % CAK sequence for pre-optimisation
	info = ss_info(A,C,K,V);    % SS information
	fres = info.fres;           % frequency resolution
	H = ss2trfun(A,C,K,fres);   % transfer function
	gc = ss_to_pwcgc(A,C,K,V);  % causal graph
end
rng_restore(mrstate);

% Display causal graph

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);

% Set 1+1 evolutionary strategy parameters

algo = '1+1 ES';
[ifac,nfac] = es_parms(esrule,m*(n-m));

irstate = rng_seed(iseed);

% Normalisation factors

if nnorm > 0
	fprintf('\n');
	L = rand_orthonormal(n,m,nnorm);
	D = zeros(nnorm,1);
	for k = 1:nnorm
		D(k) = cak2ddx(L(:,:,k),CAK);
	end
	ddpmean = mean(D);
	fprintf('dd proxy mean    = %g\n',ddpmean);
	for k = 1:nnorm
		D(k) = trfun2dd(L(:,:,k),H);
	end
	ddsmean = mean(D);
	fprintf('dd spectral mean = %g\n',ddsmean);
	for k = 1:nnorm
		D(k) = iss2dd(L(:,:,k),A,C,K);
	end
	ddnmean = mean(D);
	fprintf('dd mean          = %g\n\n',ddnmean);
end

% Initialise optimisation

iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = rand_orthonormal(n,m,nruns); % initial (orthonormalised) random linear projections

rng_restore(irstate);

if hist
	dhistp = nan(piters,nruns);
	dhists = nan(siters,nruns);
	dhistn = nan(niters,nruns);
end

orstate = rng_seed(oseed);
for k = 1:nruns

	fprintf('run %2d of %2d\n',k,nruns);

	Loptk = Lopt(:,:,k);
	sigk = sig0;

	if piters > 0 % optimisation using "proxy" DD

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddx(CAK,Loptk,piters,sigk,ifac,nfac,estol,hist);
		if hist, dhistp(:,k) = dhistk; end
		fprintf('\topt 1 : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged, fprintf('converged  '); else, fprintf('unconverged'); end
		fprintf(' in %6d iterations\n',ioptk);

	end

	if siters > 0 % optimisation using integrated spectral DD

		if sigreset, sigk = sig0; end

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dds(H,Loptk,siters,sigk,ifac,nfac,estol,hist);
		if hist, dhists(:,k) = dhistk; end
		fprintf('\topt 2 : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged, fprintf('converged  '); else, fprintf('unconverged'); end
		fprintf(' in %6d iterations\n',ioptk);

	end

	if niters > 0 % optimisation using state-space DD (DARE)

		if sigreset, sigk = sig0; end

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dd(A,C,K,Loptk,niters,sigk,ifac,nfac,estol,hist);
		if hist, dhistn(:,k) = dhistk; end
		fprintf('\topt 3 : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged, fprintf('converged  '); else, fprintf('unconverged'); end
		fprintf(' in %6d iterations\n',ioptk);

	end

	Lopt(:,:,k) = Loptk;
	dopt(k) = doptk;
	iopt(k) = ioptk;

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

% Save workspace

clear k
if hist
	wsfile = fullfile(resdir,[scriptname '_hist' rid '.mat']);
else
	clear dhistp dhistn
	wsfile = fullfile(resdir,[scriptname '_nohist' rid '.mat']);
end
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');

% Plot optimisation histories

if hist
	if nnorm > 0
		dhistp = dhistp/ddpmean;
		dhists = dhists/ddsmean;
		dhistn = dhistn/ddnmean;
	end
	gptitle = sprintf('Optimisation history (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
	gpstem = fullfile(resdir,[scriptname '_opthist' rid]);
	gp_opthist(dhistp,dhists,dhistn,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
