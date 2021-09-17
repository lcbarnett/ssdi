function [d,L,converged,sig,iters,dhist] = opt_es_gd_dds(H,L,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - L is orthonormal
% 2 - Residuals covariance matrix is identity

if isscalar(tol)
	stol = tol;
	dtol = tol;
else
	stol = tol(1);
	dtol = tol(2);
end

[n,m] = size(L);

% Calculate dynamical dependence of initial projection

d = trfun2dd(L,H);

if hist
	dhist = zeros(maxiters,3);
	i = 1; dhist(i,:) = [1 d sig];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	if hist, i = i+1; dhist(i,:) = [iters d sig]; end

	% Move down gradient and orthonormalise

	L = orthonormalise(L - sig*trfun2ddgrad(L,H)); % this is not quite correct, but probably the best we can do :-/

	% Calculate dynamical dependence of candidate projection

	dold = d;
	d = trfun2dd(L,H);

	% If dynamical dependence smaller, increase step size, else decrease

	if d < dold
		sig = ifac*sig;
	else
		sig = nfac*sig;
	end

	if hist, i = i+1; dhist(i,:) = [iters d sig]; end

	% Test convergence

	if sig < stol
		converged = 1;
		break
	end

	if d < dtol
		converged = 2;
		break
	end

end

if hist
	i = i+1; dhist(i,:) = [iters d sig];
	dhist = dhist(1:i,:);
end
