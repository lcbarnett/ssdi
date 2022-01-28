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

if ~exist('resdir',   'var'), resdir   = tempdir;     end % results directory
if ~exist('fnamep',   'var'), fnamep   = 'preopt_dd'; end % pre-optimisation filename root
if ~exist('fnameo',   'var'), fnameo   = 'opt_dd';    end % pre-optimisation filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ctol',     'var'), ctol     = 1e-6;        end % hyperplane clustering tolerance
if ~exist('hist',     'var'), hist     = true;        end % calculate optimisation history?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('sig0o',    'var'), sig0o    = 0.1;         end % optimisation (gradient descent) initial step size
if ~exist('gdlso',    'var'), gdlso    = [2,1/2];     end % gradient-descent "line search" parameters
if ~exist('gdtolo',   'var'), gdtolo   = 1e-10;       end % gradient descent convergence tolerance
if ~exist('niterso',  'var'), niterso  = 10000;       end % pre-optimisation iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';     end % Gnuplot terminal
if ~exist('gpscale',  'var'), gpscale  = [Inf,1.1];   end % Gnuplot scale
if ~exist('gpfsize',  'var'), gpfsize  = 14;          end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;           end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

% Load pre-optimisation workspace

if hist
	wsfilep = fullfile(resdir,[fnamep '_mdim_' num2str(mdim) '_H.mat']);
else
	wsfilep = fullfile(resdir,[fnamep '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** loading pre-optimisation workspace from ''%s''... ',wsfilep);
load(wsfilep);
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
	gp_opthist({histp;histo},[nitersp;niterso],[true;true],{'Pre-optimisation (GD)','Optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance: %s, m = %d',mdescript,m);
gpstem = fullfile(tempdir,'opt_iodist');
gp_iodist(gopto,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

% Save optimisation workspace

clear n m k st et tcpu gpplot gpterm gpscale gpfsize gptitle gpstem

if hist
	wsfileo = fullfile(resdir,[fnameo '_mdim_' num2str(mdim) '_H.mat']);
else
	wsfileo = fullfile(resdir,[fnameo '_mdim_' num2str(mdim) '_N.mat']);
end
fprintf('\n*** saving optimisation workspace in ''%s''... ',wsfileo);
save(wsfileo);
fprintf('done\n\n');
