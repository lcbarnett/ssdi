%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstrate calculation of Causal Emergence (CE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Supply model for your data with decorrelated and normalised residuals (see,
% e.g., sim_model.m). Model properties required for this script are:
%
% mdescript - a model description string
% H         - transfer function for DD calculation by spectral method
% G, G0     - covariance/autocovariance sequence (for sim_model, set aclmax to, say, 10000)
%
% Specify a macroscopic dimension mdim and simulation parameters, or accept
% defaults (see below).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   tempdir     );  % model directory
defvar('modname',  'sim_model' );  % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('nsamps',   100         ); % number of sample random projections
defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',   'x11'       ); % Gnuplot terminal
defvar('gpscale',  [1,1.2]     ); % Gnuplot scale
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

fprintf('\n%s: CE calculation for m = %d\n',mdescript,m);

fprintf('\nCalculating the Sigma_i ');
st = tic;
[CES,G0c] = ac2ces(G);
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Random projection

rstate = rng_seed(iseed);
L = rand_orthonormal(n,mdim,nsamps); % (orthonormalised) random linear projections
rng_restore(rstate);

npi = floor(nsamps/10); % for progress indicator

% Calculate DDs (spectral method - could alternatively use DARE method)

fprintf('\nCalculating dynamical dependences ');
st = tic;
DD = zeros(nsamps,1);
for i = 1:nsamps
	DD(i) = trfun2dd(L(:,:,i),H);
	if ~mod(i,npi), fprintf('.'); end % progress indicator
end
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Calculate CEs

fprintf('\nCalculating causal emergences ');
st = tic;
CE = zeros(nsamps,1);
for i = 1:nsamps
	CE(i) = ces2ce(L(:,:,i),CES,G0c,DD(i));
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
