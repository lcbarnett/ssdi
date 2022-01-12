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
if ~exist('ssig0',    'var'), ssig0    = 0.01;      end % spectral optimisation initial step size
if ~exist('gsig0',    'var'), gsig0    = 0.1;       end % gradient descent optimisation initial step size
if ~exist('dsig0',    'var'), dsig0    = 0.001;     end % state-space optimisation initial step size
if ~exist('esrule',   'var'), esrule   = 1/5;       end % evolution strategy step-size adaptation rule
if ~exist('estol',    'var'), estol    = 1e-8;      end % evolution strategy convergence tolerance
if ~exist('fres',     'var'), fres     = [];        end % frequency resolution (empty for automatic)
if ~exist('sitol',    'var'), sitol    = 1e-12;     end % spectral integration tolerance
if ~exist('siminp2',  'var'), siminp2  = 6;         end % spectral integration freq. res. min power of 2
if ~exist('simaxp2',  'var'), simaxp2  = 14;        end % spectral integration freq. res. max power of 2
if ~exist('nsics',    'var'), nsics    = 100;       end % number of samples for spectral integration check
if ~exist('npiters',  'var'), npiters  = 100000;    end % pre-optimisation iterations
if ~exist('nsiters',  'var'), nsiters  = 1000;      end % spectral optimisation iterations
if ~exist('ngiters',  'var'), ngiters  = 1000;      end % gradient descent optimisation iterations
if ~exist('nditers',  'var'), nditers  = 1000;      end % SS optimisation iterations
if ~exist('nruns',    'var'), nruns    = 10;        end % runs (restarts)
if ~exist('hist',     'var'), hist     = true;      end % calculate optimisation history?
if ~exist('iseed',    'var'), iseed    = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('oseed',    'var'), oseed    = 0;         end % optimisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varmod = exist('G','var');
if varmod % VAR model
	if isscalar(G), n = G; G = ones(n); else, n = size(G,1); end; % fully-connected if G scalar
	if ~exist('r',    'var'), r = n;   end % VAR model order
	if ~exist('w',    'var'), w = 1;   end % VAR coefficients decay parameter
else      % fully-connected state-space model
	if ~exist('n',    'var'), n = 9;   end % microscopic dimension
	if ~exist('r',    'var'), r = 3*n; end % hidden state dimension
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',   'var'), gvprog   = 'neato';   end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',   'var'), gvdisp   = true;      end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;
siparms    = [sitol siminp2 simaxp2];

% Generate random VAR or ISS model

rstate = rng_seed(mseed);
V0 = corr_rand(n,rcorr); % residuals covariance matrix
if varmod
	ARA0 = var_rand(G,r,rho,w);             % random VAR model
	gc = var_to_pwcgc(ARA0,V0);             % causal graph
	[ARA,V] = transform_var(ARA0,V0);       % transform model to decorrelated-residuals form
	[A,C,K] = var_to_ss(ARA);               % equivalent ISS model
	if isempty(fres)
		[fres,ierr] = var2fres(ARA,V,siparms);
		if isnan(fres) % failed!
			fprintf(2,'WARNING: Spectral integral frequency resolution estimation failed - defaulting to autocorrelation estimate');
			[fres,ierr] = var2fres(ARA,V);  % use autocorrelation-based estimate
		end
	end
	H = var2trfun(ARA,fres);                % transfer function
	CAK = ARA;                              % CAK sequence for pre-optimisation
else
	[A0,C0,K0] = iss_rand(n,r,rho);         % random ISS model
	gc = ss_to_pwcgc(A0,C0,K0,V0);          % causal graph
	[A,C,K,V] = transform_ss(A0,C0,K0,V0);  % transform model to decorrelated-residuals form
	if isempty(fres)
		[fres,ierr] = ss2fres(A,C,K,V,siparms);
		if isnan(fres) % failed!
			fprintf(2,'\nWARNING: Spectral integral frequency resolution estimation failed - defaulting to autocorrelation estimate');
			[fres,ierr] = ss2fres(A,C,K,V); % use autocorrelation-based estimate
		end
	end
	H = ss2trfun(A,C,K,fres);               % transfer function
	CAK = iss2cak(A,C,K);                   % CAK sequence for pre-optimisation
end
fprintf('\nFrequency resolution = %d (integration error = %e)\n\n',fres,ierr);
fprintf('Spectral DD accuracy check (frequency resolution = %d) ... ',fres);
derr = dds_check(A,C,K,H,m,nsics); % spectral integration check
fprintf('integration error = %e\n\n',derr);
if derr > 1e-12, fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n'); end
rng_restore(rstate);

% Display causal graph

eweight = gc/nanmax(gc(:));
gfile = fullfile(resdir,[scriptname '_pwcgc' rid]);
wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);
fprintf('\n');

% Set 1+1 evolution strategy parameters

algo = '1+1 ES';
[ifac,nfac] = es_parms(esrule,m*(n-m));

% Initialise optimisation

rstate = rng_seed(iseed);
iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = rand_orthonormal(n,m,nruns); % initial (orthonormalised) random linear projections
rng_restore(rstate);
if hist
	dhistp = cell(nruns,1);
	dhists = cell(nruns,1);
	dhistg = cell(nruns,1);
	dhistd = cell(nruns,1);
end

rstate = rng_seed(oseed);
for k = 1:nruns

	fprintf('run %2d of %2d\n',k,nruns);

	Loptk = Lopt(:,:,k);

	% Pre-optimisation using "proxy" DD

	if npiters > 0
		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddx(CAK,Loptk,npiters,psig0,ifac,nfac,estol,hist);
		if hist, dhistp{k} = dhistk; end
		fprintf('\tpre-opt : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);
	end

	% Optimisation using integrated spectral method (generally faster, potentially less accurate)

	if nsiters > 0
		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dds(H,Loptk,nsiters,ssig0,ifac,nfac,estol,hist);
		if hist, dhists{k} = dhistk; end
		fprintf('\tsopt    : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);
	end

	% Optimisation using gradient descent

	if ngiters > 0
		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_ddg(H,Loptk,ngiters,gsig0,ifac,nfac,estol,hist);
		if hist, dhistg{k} = dhistk; end
		fprintf('\tgopt    : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);
	end

	% Optimisation using state-space (DARE) method (most accurate, but may be slower)

	if nditers > 0       % use state-space (DARE) method (usually slower)
		[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dd(A,C,K,Loptk,nditers,dsig0,ifac,nfac,estol,hist);
		if hist, dhistd{k} = dhistk; end
		fprintf('\tdopt    : dopt = %.4e : sig = %.4e : ',doptk,sigk);
		if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
		fprintf(' in %6d iterations\n',ioptk);
	end

	Lopt(:,:,k) = Loptk;
	dopt(k) = doptk;
	iopt(k) = ioptk;

end
rng_restore(rstate);

% Transform Lopt back to correlated residuals form

Lopt = transform_proj(Lopt,V0); % V0 is the original residuals covariance matrix

% Sort (local) optima by dynamical dependence

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
if hist
	dhistp = dhistp(sidx);
	dhists = dhists(sidx);
	dhistg = dhistg(sidx);
	dhistd = dhistd(sidx);
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
	dhist  = {dhistp;dhists;dhistg;dhistd};
	niters = [npiters;nsiters;ngiters;nditers];
	titles = {'Pre-optimisation';'Spectral optimisation';'Gradient descent';'SS optimisation'};
	gp_opthist(dhist,niters,titles,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance (%s) : n = %d, r = %d, m = %d',algo,n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gpscale = [1.2,1.1];
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
