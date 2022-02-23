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
defvar('sig0o',    0.1         ); % optimisation (gradient descent) initial step size
defvar('gdlso',    2           ); % gradient-descent "line search" parameters
defvar('gdtolo',   1e-10       ); % gradient descent convergence tolerance
defvar('histo',    true        ); % calculate optimisation history?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x-pdf'     ); % Gnuplot terminal
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

[uidx,usiz,nrunso] = Lcluster(goptp,ctol,doptp);

% Initialise optimisations

L0o = Loptp(:,:,uidx);

% Multiple optimisation runs

st = tic;
[dopto,Lopto,convp,iopto,sopto,cputo,ohisto,et] = opt_gd_dds_mruns(H,L0o,nrunso,niterso,sig0o,gdlso,gdtolo,histo);
et = toc(st);

fprintf('\noptimal dynamical dependence =\n'); disp(dopto');
fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));
fprintf('CPU secs per run = %7.4f +- %6.4f\n\n',mean(cputo),std(cputo));

% Plot optimisation histories

if histo
	gptitle  = sprintf('Optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'opt_hist');
	gp_opthist({ohistp;ohisto},[nitersp;niterso],[true;true],true,{'Pre-optimisation (GD)','Optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gopto = gmetrics(Lopto);
gptitle = sprintf('Inter-optimum distance: %s, m = %d',mdescript,m);
gpstem = fullfile(tempdir,'opt_iodist');
gp_iodist(gopto,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

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

% !!! Now remember to convert Lopto back to non-decorrelated coordinates !!!
