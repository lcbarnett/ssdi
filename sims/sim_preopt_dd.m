%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence pre-optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model parameters (see sim_model.m). Must supply m = macroscopic state
% dimension, then run this script as, e.g.,
%
% clear; m = 4; G = tnet9d; mseed = 923933; sim_preopt_dd
%
% Now run Lcluster to get a good uidx tolerance ctol; e.g.:
%
% ctol = 1e-6; [uidx usiz] = Lcluster(goptp,ctol,doptp)
%
% followed by sim_opt_dd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('resdir',   'var'), resdir   = tempdir;   end % results directory
if ~exist('rid',      'var'), rid      = '';        end % run ID tag

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('sig0p',    'var'), sig0p    = 1;         end % pre-optimisation (gradient descent) initial step size
if ~exist('gdlsp',    'var'), gdlsp    = [2,1/2];   end % gradient-descent "line search" parameters
if ~exist('gdtolp',   'var'), gdtolp   = 1e-10;     end % gradient descent convergence tolerance
if ~exist('nitersp',  'var'), nitersp  = 10000;     end % pre-optimisation iterations
if ~exist('nrunsp',   'var'), nrunsp   = 100;       end % pre-optimisation runs (restarts)
if ~exist('hist',     'var'), hist     = true;      end % calculate optimisation history?
if ~exist('iseed',    'var'), iseed    = 0;         end % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpscale',  'var'), gpscale  = [Inf,0.6]; end % Gnuplot scale
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnamep = 'sim_preopt_dd';

% Generate random VAR or ISS model

sim_model;

% Initialise optimisation

doptp = zeros(1,nrunsp);
Loptp = zeros(n,m,nrunsp);
convp = false(1,nrunsp);
ioptp = zeros(1,nrunsp);
soptp = zeros(1,nrunsp);
cputp = zeros(1,nrunsp);
if hist
	histp = cell(nrunsp,1);
end

rstate = rng_seed(iseed);
L0p = rand_orthonormal(n,m,nrunsp); % initial (orthonormalised) random linear projections
rng_restore(rstate);

% "Proxy" DD pre-optimisation (gradient descent)

st = tic;
for k = 1:nrunsp
	fprintf('pre-opt run %4d of %4d : ',k,nrunsp);
	tcpu = cputime;
	[doptp(k),Loptp(:,:,k),convp(k),soptp(k),ioptp(k),histp{k}] = opt_gd_ddx(CAK,L0p(:,:,k),nitersp,sig0p,gdlsp,gdtolp,hist);
	cputp(k) = cputime-tcpu;
	fprintf('dd = %.4e : sig = %.4e : ',doptp(k),soptp(k));
	if convp(k) > 0, fprintf('converged(%d)',convp(k)); else, fprintf('unconverged '); end
	fprintf(' in %4d iterations : CPU secs = %6.2f\n',ioptp(k),cputp(k));
end
et = toc(st);

% Recalculate actual (i.e., not proxy) DDs

for k = 1:nrunsp
	doptp(k) = trfun2dd(Loptp(:,:,k),H);
end

% Sort everything by dynamical dependence

[doptp,sidxp] = sort(doptp);
Loptp = Loptp(:,:,sidxp);
ioptp = ioptp(sidxp);
convp = convp(sidxp);
soptp = soptp(sidxp);
cputp = cputp(sidxp);
if hist
	histp = histp(sidxp);
end
fprintf('\noptimal dynamical dependence =\n'); disp(doptp');

fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));

fprintf('CPU secs per run = %7.4f +- %6.4f\n',mean(cputp),std(cputp));
fprintf('\n');

goptp = gmetrics(Loptp);

% Plot optimisation histories

if hist
	gptitle  = sprintf('Pre-optimisation history : n = %d, r = %d, m = %d',n,r,m);
	gpstem   = fullfile(resdir,[fnamep '_opthist' rid]);
	gp_opthist({histp},nitersp,true,{'Pre-optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-preoptimum distance : n = %d, r = %d, m = %d',n,r,m);
gpstem = fullfile(resdir,[fnamep '_iodist' rid]);
gp_iodist(goptp,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

% Save workspace

clear k st et tcpu rstate gvdisp gvprog gpplot gpterm gpscale gpfsize gptitle gpstem
if hist
	wsfilep = fullfile(resdir,[fnamep '_hist' rid '.mat']);
else
	wsfilep = fullfile(resdir,[fnamep '_nohist' rid '.mat']);
end
fprintf('*** saving workspace in ''%s''... ',wsfilep);
save(wsfilep);
fprintf('done\n\n');
