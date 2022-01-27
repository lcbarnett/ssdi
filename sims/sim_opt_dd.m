%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamical dependence optimisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Must run sim_preopt_dd first!
%
% Then Lcluster to get a good uidx tolerance ctol; e.g.:
%
% ctol = 1e-6; [uidx usiz] = Lcluster(goptp,ctol,doptp)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('ctol',     'var'), ctol     = 1e-6;      end % tolerance for pre-optimisation clustering
if ~exist('sig0o',    'var'), sig0o    = 0.1;       end % optimisation (gradient descent) initial step size
if ~exist('gdlso',    'var'), gdlso    = [2,1/2];   end % gradient-descent "line search" parameters
if ~exist('gdtolo',   'var'), gdtolo   = 1e-10;     end % gradient descent convergence tolerance
if ~exist('niterso',  'var'), niterso  = 10000;     end % pre-optimisation iterations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gpterm',   'var'), gpterm   = 'x-pdf';   end % Gnuplot terminal
if ~exist('gpscale',  'var'), gpscale  = [Inf,1.1]; end % Gnuplot scale
if ~exist('gpfsize',  'var'), gpfsize  = 14;        end % Gnuplot font size
if ~exist('gpplot',   'var'), gpplot   = 2;         end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnameo = 'sim_opt_dd';

% Read in pre-optimisation workspace

fprintf('\n*** loading workspace in ''%s''... ',wsfilep);
load(wsfilep);
fprintf('done\n\n');

% Cluster pre-optimisation results (get indices of unique hyperplanes) to get initial configurations

[uidx,usiz] = Lcluster(goptp,ctol);
nrunso = length(uidx);
fprintf('Clusters = %d\n',nrunso);
for k = 1:nrunso
	fprintf('cluster at run %3d : size = %3d\n',uidx(k),usiz(k));
end
fprintf('\n');


% Initialise optimisation

dopto = zeros(1,nrunso);
Lopto = zeros(n,m,nrunso);
convo = false(1,nrunso);
iopto = zeros(1,nrunso);
sopto = zeros(1,nrunso);
cputo = zeros(1,nrunso);
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
	gptitle  = sprintf('Optimisation history : n = %d, r = %d, m = %d',n,r,m);
	gpstem   = fullfile(resdir,[fnameo '_opthist' rid]);
	gp_opthist({histp;histo},[nitersp;niterso],[true;true],{'Pre-optimisation (GD)','Optimisation (GD)'},gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);
end

% Plot inter-optima subspace distances

gptitle = sprintf('Inter-optimum distance : n = %d, r = %d, m = %d',n,r,m);
gpstem = fullfile(resdir,[fnameo '_iodist' rid]);
gp_iodist(gopto,gptitle,gpstem,gpterm,[1.2,1.1],gpfsize,gpplot);

% Save workspace

clear k st et tcpu gpplot gpterm gpscale gpfsize gptitle gpstem
if hist
	wsfileo = fullfile(resdir,[fnameo '_hist' rid '.mat']);
else
	wsfileo = fullfile(resdir,[fnameo '_nohist' rid '.mat']);
end
fprintf('*** saving workspace in ''%s''... ',wsfileo);
save(wsfileo);
fprintf('done\n\n');
