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
defvar('hist',     true        ); % calculate optimisation history?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('sig0o',    0.1         ); % optimisation (gradient descent) initial step size
defvar('gdlso',    2           ); % gradient-descent "line search" parameters
defvar('gdtolo',   1e-10       ); % gradient descent convergence tolerance
defvar('niterso',  10000       ); % pre-optimisation iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x-pdf'     ); % Gnuplot terminal
defvar('gpscale',  [Inf,1.1]   ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   2           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

% Load pre-optimisation results

if hist
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

fprintf('%s: optimisation for m = %d\n\n',mdescript,m);

fprintf('--------------------------------------------\n');
fprintf('Model                : %s\n',mdescript);
fprintf('--------------------------------------------\n');
fprintf('Dimension            : %d\n',n);
fprintf('Complexity (CAK)     : %d x %d x %d\n',size(CAK,1),size(CAK,2),size(CAK,3));
fprintf('Frequency resolution : %d\n',fres);
fprintf('--------------------------------------------\n');

% Cluster hyperplanes

[uidx,usiz,nrunso] = Lcluster(goptp,ctol,doptp);

% Initialise optimisation

dopto = zeros(1,  nrunso);
Lopto = zeros(n,m,nrunso);
convo = false(1,  nrunso);
iopto = zeros(1,  nrunso);
sopto = zeros(1,  nrunso);
cputo = zeros(1,  nrunso);
if hist
	histo = cell(nrunso,1);
end

L0o = Loptp(:,:,uidx);

% DD optimisation (gradient descent)

st = tic;
for k = 1:nrunso
	fprintf('opt run %2d of %2d : ',k,nrunso);
	tcpu = cputime;
	[dopto(k),Lopto(:,:,k),convo(k),sopto(k),iopto(k),histo{k}] = opt_gd_dds(H,L0o(:,:,k),niterso,sig0o,gdlso,gdtolo,hist);
	cputo(k) = cputime-tcpu;
	fprintf('dd = %.4e : sig = %.4e : ',dopto(k),sopto(k));
	if convo(k) > 0, fprintf('converged(%d)',convo(k)); else, fprintf('unconverged '); end
	fprintf(' in %4d iterations : CPU secs = %5.2f\n',iopto(k),cputo(k));
end
et = toc(st);

% Sort everything by dynamical dependence

[dopto,sidxo] = sort(dopto);
Lopto = Lopto(:,:,sidxo);
iopto = iopto(sidxo);
convo = convo(sidxo);
sopto = sopto(sidxo);
cputo = cputo(sidxo);
if hist
	histo = histo(sidxo);
end
fprintf('\noptimal dynamical dependence =\n'); disp(dopto');

fprintf('Simulation time = %s\n\n',datestr(seconds(et),'HH:MM:SS.FFF'));

fprintf('CPU secs per run = %7.4f +- %6.4f\n',mean(cputo),std(cputo));
fprintf('\n');

gopto = gmetrics(Lopto);

% Plot optimisation histories

if hist
	gptitle  = sprintf('Optimisation history: %s, m = %d',mdescript,m);
	gpstem   = fullfile(tempdir,'opt_hist');
	gp_opthist({histp;histo},[nitersp;niterso],[true;true],true,{'Pre-optimisation (GD)','Optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance: %s, m = %d',mdescript,m);
gpstem = fullfile(tempdir,'opt_iodist');
gp_iodist(gopto,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

% Save optimisation workspace

clear n m k st et tcpu gpplot gpterm gpscale gpfsize gptitle gpstem

if hist
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_H.mat']);
else
	optfile = fullfile(optdir,[optname '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving pre-optimisation results in ''%s''... ',optfile);
save(optfile);
fprintf('done\n\n');
