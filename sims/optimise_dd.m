%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Must run sim_preopt_dd first!
%
% Then run Lcluster until you are happy with the hyperplane clustering; e.g.:
%
% ctol = 1e-6; Lcluster(goptp,ctol,doptp);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('poptdir',  tempdir     ); % pre-optimisation directory
defvar('poptname', 'preopt_dd' ); % pre-optimisation filename root
defvar('optdir',   tempdir     ); % pre-optimisation directory
defvar('optname',  'opt_dd'    ); % pre-optimisation filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('ctol',     1e-6        ); % hyperplane clustering tolerance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('niterso',  10000       ); % pre-optimisation iterations
defvar('gdeso',    2           ); % gradient-descent ES version (1 or 2)
defvar('gdsig0o',  0.1         ); % optimisation (gradient descent) initial step size
defvar('gdlso',    2           ); % gradient-descent "line search" parameters
defvar('gdtolo',   1e-10       ); % gradient descent convergence tolerance
defvar('histo',    true        ); % calculate optimisation history?
defvar('ppo',      false       ); % parallelise multiple runs?
gdeso

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x11'       ); % Gnuplot terminal
defvar('gpscale',  [Inf,1.1]   ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   2           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

% Load pre-optimisation results

if histo
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_H.mat']);
else
	poptfile = fullfile(poptdir,[poptname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** loading pre-optimisation results from ''%s''... ',poptfile);
load(poptfile);
fprintf('done\n\n');

n = size(V0,1);
m = mdim;
fres = size(H,3)-1;

fprintf('%s: optimisation for m = %d\n',mdescript,m);

% Cluster hyperplanes

[uidx,usiz,nrunso] = Lcluster(goptp,ctol,doptp,gpterm,gpscale,gpfsize,gpplot);

% Initialise optimisations

L0o = Lp(:,:,uidx);

% Multiple optimisation runs

st = tic;
[dopto,Lo,convp,iopto,sopto,cputo,ohisto] = opt_gd_dds_mruns(H,L0o,niterso,gdeso,gdsig0o,gdlso,gdtolo,histo,ppo);
et = toc(st);

% Inverse-transform Lo back for un-decorrelated residuals

Lopto = transform_proj(Lo,V0);

fprintf('\noptimal dynamical dependence =\n'); disp(dopto');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputo),std(cputo));

% Plot optimisation histories

if histo && ~isempty(gpterm)
	gptitle  = sprintf('Optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'opt_hist');
	gp_opthist({ohistp;ohisto},[nitersp;niterso],[true;true],true,{'Pre-optimisation (GD)','Optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gopto = gmetrics(Lopto);
if ~isempty(gpterm)
	gptitle = sprintf('Inter-optimum distance: %s, m = %d',mdescript,m);
	gpstem = fullfile(tempdir,'opt_iodist');
	gp_iodist(gopto,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);
end

% Save optimisation results

clear n m st et tcpu gpplot gpterm gpscale gpfsize gptitle gpstem

if histo
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_H.mat']);
else
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving pre-optimisation results in ''%s''... ',optfile);
save(optfile);
fprintf('done\n\n');

% !!! Now remember to convert Lo back to non-decorrelated coordinates !!!
