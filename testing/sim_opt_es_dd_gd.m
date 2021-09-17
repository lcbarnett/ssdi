%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-space dynamical independence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('resdir', 'var'), resdir = tempdir; end % results directory
if ~exist('rid',    'var'), rid    = '';      end % run ID tag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',      'var'), rho      = 0.9;       end % spectral norm (< 1)
if ~exist('rcorr',    'var'), rcorr    = 0;         end % residuals correlation (actually multiinformation); 0 for no correlation
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

% Model type (and connectivity graph if VAR)

varmod = exist('G','var');
if varmod % VAR model
	if isscalar(G)
		n = G;                         % microscopic state dimension
		G = ones(n);                   % fully-connected
	else
		n = size(G,1);
	end
	if ~exist('r','var'), r = n;   end % VAR model order
	if ~exist('w','var'), w = 1;   end % VAR coefficients decay parameter
else      % fully-connected state-space model
	if ~exist('n','var'), n = 7;   end % microscopic state dimension
	if ~exist('r','var'), r = 3*n; end % hidden state dimension
end

% Generate random VAR or ISS model

rstate = rng_seed(mseed);
[V,SQRTV] = corr_rand(n,rcorr); % residuals covariance matrix (rcorr = 0 for identity matrix; i.e., no residuals correlation)
ISQRTV = SQRTV\eye(n);
if varmod
	ARA = var_rand(G,r,rho,w);
	info = var_info(ARA,V);     % VAR information
	gc = var_to_pwcgc(ARA,V);   % causal graph

	% Transform model to decorrelated-residuals form (by ISQRTV = inverse left-Cholesky factor of V)

	for k = 1:r
		ARA(:,:,k) = ISQRTV*ARA(:,:,k)*SQRTV;
	end
	[A,C,K] = var_to_ss(ARA);   % equivalent ISS model
	V = eye(n);

	CAK = ARA;                  % CAK sequence for pre-optimisation
	if isempty(fres), fres = info.fres; end % frequency resolution
	H = var2trfun(ARA,fres);    % transfer function
else
	[A,C,K] = iss_rand(n,r,rho);
	info = ss_info(A,C,K,V);    % SS information
	gc = ss_to_pwcgc(A,C,K,V);  % causal graph

	% Transform model to decorrelated-residuals form (by ISQRTV = inverse left-Cholesky factor of V)

	% A is unchanged
	C = ISQRTV*C;
	K = K*SQRTV;
	V = eye(n);

	CAK = iss2cak(A,C,K);       % CAK sequence for pre-optimisation
	if isempty(fres), fres = info.fres; end % frequency resolution
	H = ss2trfun(A,C,K,fres);   % transfer function
end
rng_restore(rstate);

% Display causal graph
%{
eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);
%}

% Set 1+1 evolutionary strategy parameters

algo = '1+1 ES';
[ifac,nfac] = es_parms(esrule,m*(n-m));

rstate = rng_seed(iseed);

% Spectral integration check

fprintf('Spectral DD integration check (frequency resolution = %d) ... ',fres);
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
	dhistp = cell(nruns,1);
	dhisto = cell(nruns,1);
end

rstate = rng_seed(oseed);
for k = 1:nruns

	fprintf('run %2d of %2d\n',k,nruns);

	Loptk = Lopt(:,:,k);

	% Pre-optimisation using "proxy" DD

	sigk = psig0;

	[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddx(CAK,Loptk,npiters,sigk,ifac,nfac,estol,hist);
	if hist, dhistp{k} = dhistk; end
	fprintf('\tpre-opt : dopt = %.4e : sig = %.4e : ',doptk,sigk);
	if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
	fprintf(' in %6d iterations\n',ioptk);

	% Optimisation

	sigk = osig0;

	if specopt % use integrated spectral method (usually faster)

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_gd_dds(H,Loptk,noiters,sigk,ifac,nfac,estol,hist);
		if hist, dhisto{k} = dhistk; end
		fprintf('\topt     : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);

	else       % use state-space (DARE) method (usually slower)

		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dd(A,C,K,Loptk,noiters,sigk,ifac,nfac,estol,hist);
		if hist, dhisto{k} = dhistk; end
		fprintf('\topt     : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);

	end

	Lopt(:,:,k) = Loptk;
	dopt(k) = doptk;
	iopt(k) = ioptk;

end
rng_restore(rstate);

% Transform Lopt back to correlated residuals form

for k = 1:nruns
	Lopt(:,:,k) = ISQRTV'*Lopt(:,:,k);
end

% Sort (local) optima by dynamical dependence

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
if hist
	dhistp = dhistp(sidx);
	dhisto = dhisto(sidx);
end
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
	gp_opthist(dhistp,npiters,dhisto,noiters,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gpscale = [1.2,1.1];
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
