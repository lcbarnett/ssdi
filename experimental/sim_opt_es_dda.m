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
if ~exist('psig0',    'var'), psig0    = 0.1;       end % pre-optimisation initial step size
if ~exist('osig0',    'var'), osig0    = 0.01;      end % optimisation initial step size
if ~exist('esrule',   'var'), esrule   = 1/5;       end % evolutionary strategy adaptation rule
if ~exist('estol',    'var'), estol    = 1e-8;      end % evolutionary strategy convergence tolerance
if ~exist('fres',     'var'), fres     = [];        end % frequency resolution (empty for automatic)
if ~exist('nsics',    'var'), nsics    = 100;       end % number of samples for spectral integration check
if ~exist('specopt',  'var'), specopt  = true;      end % use spectral DD method? (Usually faster)
if ~exist('npiters',  'var'), npiters  = 100000;    end % pre-optimisation iterations
if ~exist('noiters',  'var'), noiters  = 10000;     end % optimisation iterations
if ~exist('nruns',    'var'), nruns    = 10;        end % runs (restarts)
if ~exist('hist',     'var'), hist     = true;      end % calculate optimisation history?
if ~exist('iseed',    'var'), iseed    = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('oseed',    'var'), oseed    = 0;         end % optimisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',   'var'), gvprog   = 'neato';   end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',   'var'), gvdisp   = true;      end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;

% Generate random VAR or ISS model

rstate = rng_seed(mseed);
V = eye(n); % residuals covariance is decorrelated and normalised to unity
if varmod
	ARA = var_rand(G,r,rho,w);
	[A,C,K] = var_to_ss(ARA);
	CAK = ARA;                  % CAK sequence for pre-optimisation
	info = var_info(ARA,V);     % VAR information
	if isempty(fres), fres = info.fres; end % frequency resolution
	H = var2trfun(ARA,fres);    % transfer function
	gc = var_to_pwcgc(ARA,V);   % causal graph
else
	[A,C,K] = iss_rand(n,r,rho);
	CAK = iss2cak(A,C,K);       % CAK sequence for pre-optimisation
	info = ss_info(A,C,K,V);    % SS information
	if isempty(fres), fres = info.fres; end % frequency resolution
	H = ss2trfun(A,C,K,fres);   % transfer function
	gc = ss_to_pwcgc(A,C,K,V);  % causal graph
end
rng_restore(rstate);

% Display causal graph

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);

% Set 1+1 evolutionary strategy with evolutimg step-size parameters

algo = '1+1 ESS';
ssig = 0.5/sqrt(2*m*(n-m));

rstate = rng_seed(iseed);

% Spectral integration check

fprintf('\nSpectral DD integration check (frequency resolution = %d) ... ',fres);
L = rand_orthonormal(n,m,nsics); % (orthonormalised) random linear projections
D1 = zeros(nsics,1);
for k = 1:nsics
	D1(k) = iss2dd(L(:,:,k),A,C,K);
end
D2 = zeros(nsics,1);
for k = 1:nsics
	D2(k) = trfun2dd(L(:,:,k),H);
end
mad = maxabs(D1-D2);
fprintf('max. abs. diff = %e\n\n',mad);
if mad > 1e-12, fprintf(2,'WARNING: spectral method may be inaccurate!\n\n'); end
clear L D1 D2

% Initialise optimisation

iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = rand_orthonormal(n,m,nruns); % initial (orthonormalised) random linear projections

rng_restore(rstate);

if hist
	dhistp = nan(npiters,2,nruns);
	dhisto = nan(noiters,2,nruns);
end

rstate = rng_seed(oseed);
for k = 1:nruns

	fprintf('run %2d of %2d\n',k,nruns);

	Loptk = Lopt(:,:,k);

	% Pre-optimisation using "proxy" DD

	sigk = psig0;

	[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddxa(CAK,Loptk,npiters,sigk,ssig,estol,hist);
	if hist, dhistp(:,:,k) = dhistk; end
	fprintf('\tpre-opt : dopt = %.4e : sig = %.4e : ',doptk,sigk);
	if converged, fprintf('converged  '); else, fprintf('unconverged'); end
	fprintf(' in %6d iterations\n',ioptk);

	% Optimisation

	sigk = osig0;

	if specopt % use integrated spectral method (usually faster)

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddsa(H,Loptk,noiters,sigk,ssig,estol,hist);
		if hist, dhisto(:,:,k) = dhistk; end
		fprintf('\topt     : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged, fprintf('converged  '); else, fprintf('unconverged'); end
		fprintf(' in %6d iterations\n',ioptk);

	else       % use state-space (DARE) method (usually slower)

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dda(A,C,K,Loptk,noiters,sigk,ssig,estol,hist);
		if hist, dhisto(:,:,k) = dhistk; end
		fprintf('\topt     : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged, fprintf('converged  '); else, fprintf('unconverged'); end
		fprintf(' in %6d iterations\n',ioptk);

	end

	Lopt(:,:,k) = Loptk;
	dopt(k) = doptk;
	iopt(k) = ioptk;

end
rng_restore(rstate);

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

clear k doptk Loptk sigk ioptk dhistk converged
if hist
	wsfile = fullfile(resdir,[scriptname '_hist' rid '.mat']);
else
	wsfile = fullfile(resdir,[scriptname '_nohist' rid '.mat']);
end
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');

% Plot optimisation histories

if hist
	gptitle = sprintf('Optimisation history (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
	gpstem = fullfile(resdir,[scriptname '_opthist' rid]);
	gpscale = [Inf,1.5];
	gp_opthist(dhistp,dhisto,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gpscale = [1.2,1.1];
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
