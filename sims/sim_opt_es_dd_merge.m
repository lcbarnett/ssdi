% Specify mergedir, mergeroot, mergerid
%
% probably
%
% mergeroot = 'sim_opt_es_dd_nohist';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',  'var'), gvprog  = 'dot';   end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',  'var'), gvdisp  = true;    end % GraphViz display? (else just generate graph files)
if ~exist('gpterm',  'var'), gpterm  = 'x-pdf'; end % Gnuplot terminal
if ~exist('gpscale', 'var'), gpscale = 1.2;     end % Gnuplot scale factor(s)
if ~exist('gpfsize', 'var'), gpfsize = 14;      end % Gnuplot font size
if ~exist('gpplot',  'var'), gpplot  = 2;       end % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mm = 1;
matfile = fullfile(mergedir,[mergeroot mergerid '_' num2str(mm) '.mat']);
load(matfile);

n1 = n-1;

tmp_iopt = zeros(nruns,n1);
tmp_dopt = zeros(nruns,n1);
tmp_Lopt = cell(n1,1);

tmp_iopt(:,mm) = iopt;
tmp_dopt(:,mm) = dopt;
tmp_Lopt{mm} = Lopt;

for mm = 2:n1
	matfile = fullfile(mergedir,[mergeroot mergerid '_' num2str(mm) '.mat']);
	fprintf('merging ''%s''\n',matfile);
	load(matfile);
	tmp_iopt(:,mm) = iopt;
	tmp_dopt(:,mm) = dopt;
	tmp_Lopt{mm} = Lopt;
end

iopt = tmp_iopt;
dopt = tmp_dopt;
Lopt = tmp_Lopt;

resdir     = mergedir;
rid        = mergerid;
scriptname = mfilename;

clear mergedir mergeroot mergerid tmp_iopt tmp_dopt tmp_Lopt mm

% fprintf('\noptimal dynamical dependence =\n\n');
% disp(num2str(1:n1,'    %8d'));
% disp(num2str(dopt,'    %8.6f'));

eweight = gc/nanmax(gc(:));

nweight = zeros(n,n1,nruns);
for m = 1:n1
	for k = 1:nruns
		nweight(:,m,k) = 1-gmetricsx(Lopt{m}(:,:,k));
	end
end

Loptd = zeros(nruns,nruns,n1);
for m = 1:n1
	Loptd(:,:,m) = gmetrics(Lopt{m});
end

%%{
Loptx = cell(n1,1);
for m = 1:n1
	fprintf('m = %d\n',m);
	Loptx{m} = zeros(nchoosek(n,m),nruns);
	for k = 1:nruns
		Loptx{m}(:,k) = 1-gmetricsxx(Lopt{m}(:,:,k));
	end
end
%%}

clear m k

wsfile = fullfile(resdir,[scriptname rid '.mat']);
fprintf('*** saving workspace in ''%s''... ',wsfile);
save(wsfile);
fprintf('done\n');

% Plot (local) optimum dynamical dependencies at all scales

gptitle = sprintf('Local optima (%s) : n = %d, r = %d',algo,n,r);
gpstem = fullfile(resdir,[scriptname '_localopt' rid]);
gp_localopt(dopt,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

gptitle = sprintf('Inter-subspace distances (%s) : n = %d, r = %d',algo,n,r);
gpstem = fullfile(resdir,[scriptname '_iodist' rid]);
gp_iodist_all(Loptd,gptitle,gpstem,gpterm,gpscale,gpfsize,gpplot);

% to display Plucker distances:

gptitle = sprintf('Plucker distances (%s) : n = %d, r = %d',algo,n,r);
gpstem = fullfile(resdir,[scriptname '_plucker' rid]);

% m = 4; k = 3; gp_plucker(Loptx{m}(:,k),n,m,gptitle,fullfile(resdir,sprintf('%s_plucker%s_scale%02d_run%03d',scriptname,rid,m,k)),gpterm,gpscale,gpfsize,gpplot);


% To display axis projection-weighted graph, e.g.:

% m = 4; k = 3; wgraph2dot(nweight(:,m,k),eweight,fullfile(resdir,sprintf('%s_graph%s_scale%02d_run%03d',scriptname,rid,m,k)),[],gvprog);
