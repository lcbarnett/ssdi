%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('n','var'), n = 7;   end % microscopic state dimension
if ~exist('r','var'), r = 3*n; end % hidden state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',     'var'), rho     = 0.9;       end % spectral norm (< 1)
if ~exist('fres',    'var'), fres    = [];        end % frequency resolution (empty for auto)
if ~exist('mseed',   'var'), mseed   = 0;         end % model random seed (0 to use current rng state)
if ~exist('iseed',   'var'), iseed   = 0;         end % initialisation random seed (0 to use current rng state)
if ~exist('nsamps',  'var'), nsamps  = 1000;      end % number of samples

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

% Random projections

rstate = rng_seed(iseed);
L = randn(n,m,nsamps);
rng_restore(rstate);
for k = 1:nsamps
	L(:,:,k) = orthonormalise(L(:,:,k));
end

D1 = zeros(nsamps,1);
ptic
for k = 1:nsamps
	D1(k) = iss2dd(L(:,:,k),A,C,K);
end
ptoc

D2 = zeros(nsamps,1);
ptic
for k = 1:nsamps
	D2(k) = trapz(trfun2sdd(L(:,:,k),H))/fres;
end
ptoc

fprintf('\nmax absolute difference = %e\n',maxabs(D1-D2));
