%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Must supply m = macroscopic state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('psig0',    'var'), psig0    = 1;         end % pre-optimisation (gradient descent) initial step size
if ~exist('ssig0',    'var'), ssig0    = 0.1;       end % dd optimisation (gradient descent) initial step size
if ~exist('dsig0',    'var'), dsig0    = 0.0001;    end % dd optimisation (evolution strategy) initial step size
if ~exist('gdrule',   'var'), gdrule   = [2,1/2];   end % gradient-descent "line search" parameters
if ~exist('gdtol',    'var'), gdtol    = 1e-10;     end % gradient descent convergence tolerance
if ~exist('esrule',   'var'), esrule   = 1/5;       end % evolution strategy step-size adaptation rule
if ~exist('estol',    'var'), estol    = 1e-8;      end % evolution strategy convergence tolerance
if ~exist('npiters',  'var'), npiters  = 10000;     end % pre-optimisation (gradient descent) iterations
if ~exist('nsiters',  'var'), nsiters  = 1000;      end % spectral method optimisation (gradient descent) iterations
if ~exist('nditers',  'var'), nditers  = 100;       end % state-space method optimisation (evolution strategy) iterations
if ~exist('nruns',    'var'), nruns    = 10;        end % runs (restarts)
if ~exist('hist',     'var'), hist     = true;      end % calculate optimisation history?
if ~exist('iseed',    'var'), iseed    = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('oseed',    'var'), oseed    = 0;         end % optimisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('resdir',   'var'), resdir   = tempdir;   end % results directory
if ~exist('rid',      'var'), rid      = '';        end % run ID tag
if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scriptname = mfilename;

% Generate random VAR or ISS model

sim_model;

% Set gradient descent strategy parameters

ifgd = gdrule(1);
nfgd = gdrule(2);

% Set 1+1 evolution strategy parameters

[ifes,nfes] = es_parms(esrule,m*(n-m));

% Initialise optimisation

rstate = rng_seed(iseed);
iopt = zeros(1,nruns);
dopt = zeros(1,nruns);
Lopt = rand_orthonormal(n,m,nruns); % initial (orthonormalised) random linear projections
rng_restore(rstate);
if hist
	dhistp = cell(nruns,1);
	dhists = cell(nruns,1);
	dhistd = cell(nruns,1);
end

rstate = rng_seed(oseed);
st = tic;
for k = 1:nruns

	fprintf('run %2d of %2d\n',k,nruns);

	Loptk = Lopt(:,:,k);

	% "Proxy" DD pre-optimisation (gradient descent)

	tcpu = cputime;
	[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_gd_ddx(CAK,Loptk,npiters,psig0,ifgd,nfgd,gdtol,hist);
	secs = cputime-tcpu;
	if hist, dhistp{k} = dhistk; end
	fprintf('\tpopt : dd = %.4e : sig = %.4e : ',doptk,sigk);
	if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
	fprintf(' in %4d iterations (%.2f secs)\n',ioptk,secs);

	% DD optimisation (gradient descent) using spectral integration method

	tcpu = cputime;
	[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_gd_dds(H,Loptk,nsiters,ssig0,ifgd,nfgd,gdtol,hist);
	secs = cputime-tcpu;
	if hist, dhists{k} = dhistk; end
	fprintf('\tsopt : dd = %.4e : sig = %.4e : ',doptk,sigk);
	if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
	fprintf(' in %4d iterations (%.2f secs)\n',ioptk,secs);

	% DD optimisation (evolutionary strategy) using state-space (DARE) method (most accurate, but may be slower)

	tcpu = cputime;
	[doptk,Loptk,converged,sigk,ioptk,dhistk] = opt_es_dd(A,C,K,Loptk,nditers,dsig0,ifes,nfes,estol,hist);
	secs = cputime-tcpu;
	if hist, dhistd{k} = dhistk; end
	fprintf('\tdopt : dd = %.4e : sig = %.4e : ',doptk,sigk);
	if converged > 0, fprintf('converged(%d)',converged); else, fprintf('unconverged '); end
	fprintf(' in %4d iterations (%.2f secs)\n',ioptk,secs);

	Lopt(:,:,k) = Loptk;
	dopt(k) = doptk;
	iopt(k) = ioptk;

end
et = toc(st);
rng_restore(rstate);

fprintf('\ntotal time = %s\n',datestr(seconds(et),'HH:MM:SS.FFF'));

% Transform Lopt back to correlated residuals form

Lopt = transform_proj(Lopt,V0); % V0 is the original residuals covariance matrix

% Sort (local) optima by dynamical dependence

[dopt,sidx] = sort(dopt);
iopt = iopt(sidx);
Lopt = Lopt(:,:,sidx);
if hist
	dhistp = dhistp(sidx);
	dhists = dhists(sidx);
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
	gptitle = sprintf('Optimisation history : n = %d, r = %d, m = %d',n,r,m);
	gpstem   = fullfile(resdir,[scriptname '_opthist' rid]);
	gpscale  = [Inf,1.5];
	dhist    = {dhistp;dhists;dhistd};
	niters   = [npiters;nsiters;nditers];
	havegrad = [true,true,false];
	titles   = {'Pre-optimisation (GD)';'Spectral optimisation (GD)';'SS optimisation (ES)'};
	gp_opthist(dhist,niters,havegrad,titles,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance : n = %d, r = %d, m = %d',n,r,m);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gpscale = [1.2,1.1];
gp_iodist(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% To display axis projection-weighted graph, e.g.:
%
% k = 3; wgraph2dot(nweight(:,k),eweight,fullfile(resdir,sprintf('graph%s_run%03d',rid,k)),[],gvprog);
