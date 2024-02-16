%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demonstrate calculation of coinformation (CI) and Causal Emergence (CE)
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

defvar('nsamps',   100         ); % number of sample random projections
defvar('iseed',    0           ); % initialisation random seed (0 to use current rng state)
defvar('precomp',  true        ); % precompute quantities which don't depend on projection (speeds up computation)?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gpterm',  'x11'        ); % Gnuplot terminal
defvar('gpscale',  [1.5,0.65]  ); % Gnuplot scale
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

% Random projections

rstate = rng_seed(iseed);
L = rand_orthonormal(n,mdim,nsamps); % random linear projections (orthonormalised)
rng_restore(rstate);

% Optionally precompute quantities which don't depend on projection (speeds up computation considerably)

if precomp
	[G0,P0] = iss2ce_precomp(A0,C0,K0,V0);
else
	G0 = [];
	P0 = [];
end

% Calculate CI and DD

CI = zeros(nsamps,1);
DD = zeros(nsamps,1);
npi = floor(nsamps/10); % for progress indicator

fprintf('Calculating co-information and dynamical dependence ');
st = tic;
for i = 1:nsamps
	[CI(i),DD(i)] = iss2ce(L(:,:,i),A0,C0,K0,V0,G0,P0);
	if ~mod(i,npi), fprintf('.'); end % progress indicator
end
et = toc(st);
fprintf(' completed in %g seconds\n\n',et);

% Calculate CE

CE = CI - DD;

% Plot

gpname = 'CI_vs_DD';
gpstem = fullfile(tempdir,'CI_vs_DD');
gp_write(gpstem,[DD,CI,CE]);
gp = gp_open(gpstem,gpterm,gpscale,gpfsize);
fprintf(gp,'datfile = "%s.dat"\n\n',gpname);
fprintf(gp,'unset key\n');
fprintf(gp,'set grid\n');
fprintf(gp,'set xlabel "DD"\n');
fprintf(gp,'set xr [0:]\n');
fprintf(gp,'set multiplot title "%s: macrovariable dimension = %d\\n" layout 2,1\n',mdescript,mdim);
fprintf(gp,'set title "Co-information vs Dynamical Dependence"\n');
fprintf(gp,'set ylabel "CI" norot\n');
fprintf(gp,'plot datfile u 1:2 w points pt 7 ps 1 not\n');
fprintf(gp,'set title "Causal Emergence vs Dynamical Dependence"\n');
fprintf(gp,'set ylabel "CE" norot\n');
fprintf(gp,'plot datfile u 1:3 w points pt 7 ps 1 not\n');
fprintf(gp,'unset multiplot\n');
gp_close(gp,gpstem,gpterm,gpplot);
