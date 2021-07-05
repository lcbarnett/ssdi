function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_ddx(CAK,Lopt,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

[n,m] = size(Lopt);

% Calculate proxy dynamical dependence of initial projection

dopt = cak2ddx(Lopt,CAK);

if hist
	dhist = nan(maxiters,2);
	dhist(1) = dopt;
	dhist(2) = sig;
else
	dhist = [];
end

% Optimise

converged = false;
for iters = 2:maxiters

	% "Mutate" projection and orthonormalise

	Ltry = orthonormalise(Lopt + sig*randn(n,m));

	% Calculate proxy dynamical dependence of mutated projection

	dtry = cak2ddx(Ltry,CAK);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sig = ifac*sig;
	else
		sig = nfac*sig;
	end

	if hist
		dhist(iters,1) = dopt;
		dhist(iters,2) = sig;
	end

	% Test convergence

	if sig < tol || dopt < tol
		converged = true;
		break
	end

end
