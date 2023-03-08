
% Default parameters (override on command line - see 'defvar.h')

defvar('n',      17      ); % number of (microscopic) observables
defvar('r',      20      ); % state-space dimension
defvar('rho',    0.9     ); % AR spectral norm
defvar('g',      1       ); % innovations multi-information
defvar('mseed',  []      ); % ISS model seed (empty for no seeding)
defvar('almax',  10000   ); % maximum autocovariance lags
defvar('m',      7       ); % macroscopic dimension
defvar('lseed',  []      ); % random projection seed (empty for no seeding)

% Random ISS model

if ~isempty(mseed), rstate = rng(mseed); end
[A,C,K,rhob] = iss_rand(n,r,rho);
V = corr_rand(n,g);
if ~isempty(mseed), rng(rstate); end
fprintf('\nMA spectral norm = %g\n',rhob);

% Autocovariance sequence

[G,p] = ss_to_autocov(A,C,K,V,almax);
G0 = G(:,:,1);

fprintf('\nAutocovariance lags = %d\n',p);

fprintf('\nCalculating the Sigma_i ');
st = tic;
[CES,G0c] = ac2ces(G);
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Random projection

if ~isempty(lseed), rstate = rng(lseed); end
L = rand_orthonormal(n,m); % (orthonormalised) random linear projections
if ~isempty(lseed), rng(lseed); end

% Calculate DD (DARE method)

DD = iss2dd(L,A,C,K);

% Calculate CE (using DARE method for DD)

CE = ces2ce(L,CES,G0c,DD);

fprintf('\nDD = % g\n',  DD);
fprintf(  'CE = % g\n\n',CE);
