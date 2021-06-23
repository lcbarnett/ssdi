%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% must supply m = macroscopic state dimension

if ~exist('n','var'), n = 7;   end % microscopic state dimension
if ~exist('r','var'), r = 3*n; end % hidden state dimension

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('rho',     'var'), rho     = 0.9;       end % spectral norm (< 1)
if ~exist('mseed',   'var'), mseed   = 0;         end % model random seed (0 to use current rng state)
if ~exist('minp2',   'var'), minp2   = 4;         end % minimum power of 2
if ~exist('maxp2',   'var'), maxp2   = 16;        end % maximum power of 2
if ~exist('gpterm',  'var'), gpterm  = 'x-pdf';   end % Gnuplot terminal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random ISS model

rstate = rng_seed(mseed);
[A,C,K] = iss_rand(n,r,rho);
rng_restore(rstate);
V = eye(n);

info = ss_info(A,C,K,V);

LDV = logdet(V);

% Frequency resolution

fres = 2^nextpow2(info.acdec); % reasonable value
LDS1 = trapz(ss2ldcpsd(A,C,K,V,fres))/fres - LDV;
abserr = abs(LDS1);
fprintf('Frequency resolution estimate 1 = %d (abs. err. = %e)\n',fres,abserr);

for i = minp2:maxp2
	fres = 2^i;
	frs(i) = fres;
	LDS2(i) = trapz(ss2ldcpsd(A,C,K,V,fres))/fres - LDV;
	abserr = abs(LDS2(i));
	if abserr < eps, break; end
end
fprintf('Frequency resolution estimate 2 = %d (abs. err. = %e)\n\n',fres,abserr);

gp_qplot(frs,LDS2,[],'set title "Optimal frequency resolution estimation"\nset logs x 2\nset format x "$2^{%%L}$"\nset grid\nset xlabel "frequency resolution"\nset ylabel "error" rot\nunset key',gpterm);

function LDS = ss2ldcpsd(A,C,K,V,fres)

	[n,r,L] = ss_parms(A,C,K,V);
	h  = fres+1;
	In = eye(n);
	Ir = eye(r);
	LDS  = zeros(h,1);
	w  = exp(1i*pi*((0:fres)/fres));
	for k = 1:h % over [0,pi]
		HLk = (In + C*((w(k)*Ir-A)\K))*L;
		LDS(k) = logdet(HLk*HLk');
	end

end
