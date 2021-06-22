%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('n','var'), n = 7;   end % microscopic state dimension
if ~exist('r','var'), r = 3*n; end % hidden state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',     'var'), rho     = 0.9;       end % spectral norm (< 1)
if ~exist('fres',    'var'), fres    = [];        end % frequency resolution (empty for auto)
if ~exist('mseed',   'var'), mseed   = 0;         end % model random seed (0 to use current rng state)
if ~exist('iseed',   'var'), iseed   = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('gpterm',  'var'), gpterm  = 'x-pdf';   end % Gnuplot terminal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random ISS model

rstate = rng_seed(mseed);
[A,C,K] = iss_rand(n,r,rho);
rng_restore(rstate);

info = ss_info(A,C,K);

% Frequency resolution

if isempty(fres)
	fres = 2^nextpow2(info.acdec); % reasonable value
end
fprintf('Using frequency resolution %d\n\n',fres);

% Calculate transfer function

H = ss2trfun(A,C,K,fres);

% Random projection

rstate = rng_seed(iseed);
L = randn(n,m);
rng_restore(rstate);

L = orthonormalise(randn(n,m));

% Time-domain DD

D = iss2dd(L,A,C,K);

% Frequency-domain (spectral) DD

d = trfun2sdd(L,H);

% Plot

f = linspace(0,pi,fres+1)';
gp_qplot(f,d,[],'unset key\nset title "Spectral dynamical dependence\nset xlabel "Angular frequency"\nset ylabel "Dynamical dependence" rot\nset xr [0:pi]',gpterm);

% Check that spectral DD integrates to time-domain dd

D1 = bandlimit(d); % i.e., across full bandwidth [0,pi]

fprintf('\nIntegration check: absolute difference = %e\n',maxabs(D-D1));
