
% Default parameters (override on command line - see 'defvar.h')

defvar('n',      17      ); % number of (microscopic) observables
defvar('p',      11      ); % VAR lags
defvar('rho',    0.9     ); % VAR spectral norm
defvar('w',      1       ); % VAR decay
defvar('g',      1       ); % residuals multi-information
defvar('mseed',  []      ); % ISS model seed (empty for no seeding)
defvar('fres',   []      ); % frequency resolution (empty for automatic)
defvar('almax',  10000   ); % maximum autocovariance lags
defvar('m',      7       ); % macroscopic dimension
defvar('lseed',  []      ); % random projection seed (empty for no seeding)

% Random ISS model

if ~isempty(mseed), rstate = rng(mseed); end
A = var_rand(n,p,rho,w);
V = corr_rand(n,g);
if ~isempty(mseed), rng(rstate); end
fprintf('\nAR spectral norm = %g\n',rho);
if isempty(fres)
	[fres,ierr] = var2fres(A,V);
	fprintf('\nFrequency resolution = %d (integration error = %e)\n',fres,ierr);
end
H = var2trfun(A,fres); % transfer function

% Autocovariance sequence

[Gamma,p] = var_to_autocov(A,V,almax); % note: Gamma(:,:,1) is Gamma_0
fprintf('\nAutocovariance lags = %d\n',p);

fprintf('\nCalculating the Sigma_i ');
st = tic;
CES = ac2ces(Gamma);
et = toc(st);
fprintf(' completed in %g seconds\n',et);

% Random projection

if ~isempty(lseed), rstate = rng(lseed); end
L = rand_orthonormal(n,m); % (orthonormalised) random linear projections
if ~isempty(lseed), rng(lseed); end

% Calculate CE (using specral method for DD)

CEH = zeros(n,1);
for i = 1:n
	CEH(i) = logdet(L'*CES(:,:,i)*L);
end
DD = trfun2dd(L,H);
LG0 = L'*chol(Gamma(:,:,1),'lower');
CE = -(n-1)*logdet(LG0*LG0') - DD + sum(CEH);

fprintf('\nDD = %g\n',  DD);
fprintf(  'CE = %g\n\n',CE);
