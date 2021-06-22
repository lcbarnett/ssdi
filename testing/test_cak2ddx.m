%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('n','var'), n = 7;   end % microscopic state dimension
if ~exist('r','var'), r = 3*n; end % hidden state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',     'var'), rho     = 0.9;       end % spectral norm (< 1)
if ~exist('mseed',   'var'), mseed   = 0;         end % model random seed (0 to use current rng state)
if ~exist('iseed',   'var'), iseed   = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('nsamps',  'var'), nsamps  = 1000;      end % number of samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random ISS model

rstate = rng_seed(mseed);
[A,C,K] = iss_rand(n,r,rho);
rng_restore(rstate);

% Calculate CAK sequence for pre-optimisation

CAK = iss2cak(A,C,K);

%{
L = randn(n,m);
[L,M] = orthonormalise(L);
D1 = cak2ddx(L,M,CAK);
D2 = cak2ddx_alt(L,CAK);
%}

% Benchmark cak2ddx_alt vs cak2ddx

rstate = rng_seed(iseed);
L = randn(n,m,nsamps);
rng_restore(rstate);

fprintf('\n');

D1 = zeros(nsamps,1);
ptic
for k = 1:nsamps
	Lk = orthonormalise(L(:,:,k));
	D1(k) = cak2ddx(Lk,CAK);
end
ptoc

D2 = zeros(nsamps,1);
ptic
for k = 1:nsamps
	[Lk,Mk] = orthonormalise(L(:,:,k));
	D2(k) = cak2ddx_alt(Lk,Mk,CAK);
end
ptoc

abs(mean(D1)- mean(D2))
