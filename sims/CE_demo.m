%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstrate calculation of Causal Emergence (CE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data (see, e.g., sim_model.m). Do not transform model to
% decorrelated residuals form!
%
% Specify a macroscopic dimension mdim and simulation parameters, or accept
% defaults (see below).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   tempdir     );  % model directory
defvar('modname', 'sim_model'  );  % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('aclmax',   10000       ); % maximum autocovariance lags
defvar('nsamps',   100         ); % number of sample random projections
defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',  'x11'        ); % Gnuplot terminal
defvar('gpscale',  [1,1.2]     ); % Gnuplot scale
defvar('gpfsize',  14          ); % Gnuplot font size
defvar('gpplot',   2           ); % Gnuplot display? (0 - generate command files, 1 - generate image files, 2 - plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

assert(exist('mdim','var'),'Must supply macro dimension ''mdim''');

m = mdim;

% Load model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('\n*** loading model from ''%s''... ',modfile);
load(modfile);
fprintf('done\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: all calculations in UNTRANSFORMED coordinates!!! %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate transfer function

if varmod
	H = var2trfun(ARA0,fres);
else
	H = ss2trfun(A0,C0,K0,fres);
end

if varmod
	[G,acl] = var_to_autocov(ARA0,V0,aclmax);
else
	[G,acl] = ss_to_autocov(A0,C0,K0,V0,aclmax);
end

fprintf('\n%s: CE calculation (fres = %d, aclags = %d) for m = %d\n',mdescript,fres,acl,m);

fprintf('\nCalculating the Sigma_i ');
st = tic;
[CESRC,CLC] = ac2ces(G);
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Random projections

rstate = rng_seed(iseed);
L = rand_orthonormal(n,mdim,nsamps); % (orthonormalised) untransformed random linear projections
rng_restore(rstate);

% Calculate CEs

CE = zeros(nsamps,1);
DD = zeros(nsamps,1);
VRC = chol(V0);
npi = floor(nsamps/10); % for progress indicator

fprintf('\nCalculating causal emergence and dynamical dependence ');
st = tic;
for i = 1:nsamps
	[CE(i),DD(i)] = ces2ce(L(:,:,i),H,VRC,CESRC,CLC);
	if ~mod(i,npi), fprintf('.'); end % progress indicator
end
et = toc(st);
fprintf(' completed in %g seconds\n\n',et);

% Plot

gpstem   = fullfile(tempdir,'CE_vs_DD');
[~,gpname] = fileparts([gpstem '.xxx']); % hack to get fileparts to behave itself
gp_write(gpstem,[DD CE]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'set title "%s: CE vs DD for m = %d"\n',mdescript,mdim);
fprintf(gp,'unset key\n');
fprintf(gp,'set grid\n');
fprintf(gp,'set xlabel "DD"\n');
fprintf(gp,'set ylabel "CE" norot\n');
fprintf(gp,'plot datfile u 1:2 w points pt 7 ps 1 not\n');
gp_close(gp,gpstem,gpterm,gpplot);
