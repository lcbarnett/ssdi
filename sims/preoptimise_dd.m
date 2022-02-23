%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence pre-optimisation (via proxy DD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data with decorrelated and normalised residuals (see,
% e.g., sim_model.m). Model properties required for this script are:
%
% mdescript - a model description string
% V0        - original (i.e., pre-decorrelation) residuals covariance matrix
% CAK       - CAK sequence for pre-optimisation by proxy DD
% H         - transfer function for DD calculation by spectral method
%
% Specify a macroscopic dimension mdim and simulation parameters, or accept
% defaults (see below). After running this script, run Lcluster until you are
% happy with the hyperplane clustering; e.g.:
%
% ctol = 1e-6; Lcluster(goptp,ctol,doptp);
%
% Then run optimise_dd.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   tempdir     );  % model directory
defvar('modname',  'sim_model' );  % model filename root
defvar('poptdir',  tempdir     );  % pre-optimisation directory
defvar('poptname', 'preopt_dd' );  % pre-optimisation filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('nrunsp',   100         ); % pre-optimisation runs (restarts)
defvar('hist',     true        ); % calculate optimisation history?
defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('sig0p',    1           ); % pre-optimisation (gradient descent) initial step size
defvar('gdlsp',    2           ); % gradient-descent "line search" parameters
defvar('gdtolp',   1e-10       ); % gradient descent convergence tolerance
defvar('nitersp',  10000       ); % pre-optimisation iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x-pdf'     ); % Gnuplot terminal
defvar('gpscale',  [Inf,0.6]   ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   2           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

% Load model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('\n*** loading model from ''%s''... ',modfile);
load(modfile);
fprintf('done\n\n');

n = size(V0,1);
m = mdim;
fres = size(H,3)-1;

fprintf('%s: pre-optimisation for m = %d\n\n',mdescript,m);

fprintf('--------------------------------------------\n');
fprintf('Model                : %s\n',mdescript);
fprintf('--------------------------------------------\n');
fprintf('Dimension            : %d\n',n);
fprintf('Complexity (CAK)     : %d x %d x %d\n',size(CAK,1),size(CAK,2),size(CAK,3));
fprintf('Frequency resolution : %d\n',fres);
fprintf('--------------------------------------------\n\n');

% Initialise optimisation

doptp = zeros(1,  nrunsp);
Loptp = zeros(n,m,nrunsp);
convp = false(1,  nrunsp);
ioptp = zeros(1,  nrunsp);
soptp = zeros(1,  nrunsp);
cputp = zeros(1,  nrunsp);
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
	gptitle  = sprintf('Pre-optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'preopt_hist');
	gp_opthist({histp},nitersp,true,true,{'Pre-optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-preoptimum distance: %s, m = %d',mdescript,m);
gpstem = fullfile(tempdir,'preopt_iodist');
gp_iodist(goptp,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

% Save pre-optimisation workspace

clear n m k st et tcpu rstate gpplot gpterm gpscale gpfsize gptitle gpstem
if hist
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_H.mat']);
else
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving pre-optimisation results in ''%s''... ',poptfile);
save(poptfile);
fprintf('done\n\n');
