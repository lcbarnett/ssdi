%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grassmannian parametrisation: Stiefel vs involution (Lai-Lim-Ye, "Q-matrix")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data with decorrelated and normalised residuals (see,
% e.g., sim_model.m). Model properties required for this script are:
%
% mdescript - a model description string
% V0        - original (i.e., pre-decorrelation) residuals covariance matrix
% H         - transfer function for DD calculation by spectral method
%
% Specify a macroscopic dimension mdim and simulation parameters, or accept
% defaults (see below). After running this script, run Lcluster until you are
% happy with the hyperplane clustering; e.g.:
%
% ctol = 1e-6; Lcluster(gopts,ctol,dopts);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   tempdir     );  % model directory
defvar('modname',  'sim_model' );  % model filename root
defvar('optdir',   tempdir     );  % optimisation directory
defvar('optname',  'preopt_dd' );  % optimisation filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('nruns',    100         ); % optimisation runs (restarts)
defvar('niters',   10000       ); % optimisation iterations
defvar('hist',     true        ); % calculate optimisation history?
defvar('sig0s',    1           ); % Stiefel optimisation (gradient descent) initial step size
defvar('gdlss',    2           ); % Stiefel gradient-descent "line search" parameters
defvar('gdtols',   1e-10       ); % Stiefel gradient descent convergence tolerance
defvar('sig0i',    1           ); % LLY optimisation (gradient descent) initial step size
defvar('gdlsi',    2           ); % LLY gradient-descent "line search" parameters
defvar('gdtoli',   1e-10       ); % LLY gradient descent convergence tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x11'       ); % Gnuplot terminal
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

fprintf('%s: optimisation for m = %d\n\n',mdescript,m);

% Initialise

rstate = rng_seed(iseed);
L0 = rand_orthonormal(n,m,nruns); % initial (orthonormalised) random linear projections
rng_restore(rstate);

% Multiple Stiefel optimisation runs

st = tic;
[dopts,Ls,convs,iopts,sopts,cputs,ohists] = opt_gd_dds_mruns(H,L0,niters,sig0s,gdlss,gdtols,hist);
et = toc(st);

% Inverse-transform Lopts back for un-decorrelated residuals

Lopts = transform_proj(Ls,V0);

fprintf('\noptimal dynamical dependence =\n'); disp(dopts');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputs),std(cputs));

% Plot optimisation histories

if hist && ~isempty(gpterm)
	gptitle  = sprintf('Stiefel optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'sopt_hist');
	gp_opthist({ohists},niters,true,true,{'Stiefel optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gopts = gmetrics(Lopts);
if ~isempty(gpterm)
	gptitle = sprintf('Inter-Stiefel optimum distance: %s, m = %d',mdescript,m);
	gpstem = fullfile(tempdir,'sopt_iodist');
	gp_iodist(gopts,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);
end

% Save pre-optimisation results

clear n m st et tcpu rstate gpplot gpterm gpscale gpfsize gptitle gpstem
if hist
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_H.mat']);
else
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving optimisation results in ''%s''... ',optfile);
save(optfile);
fprintf('done\n\n');
