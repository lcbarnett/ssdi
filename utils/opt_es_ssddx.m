function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_ssddx(CAK,P0,maxiters,sig,ifac,nfac,tol,hist)

[n,m] = size(P0);

% Orthonormalise initial projection

[Lopt,Mopt] = orthonormalise(P0);

% Calculate proxy dynamical dependence of initial projection

dopt = cak2ddx(Lopt,Mopt,CAK);

if hist
	dhist = nan(maxiters,1);
	dhist(1) = dopt;
else
	dhist = [];
end

% Optimise

converged = false;
for iters = 2:maxiters

	% "Mutate" projection and orthonormalise

	[Ltry,Mtry] = orthonormalise(Lopt + sig*randn(n,m));

	% Calculate proxy dynamical dependence of mutated projection

	dtry = cak2ddx(Ltry,Mtry,CAK);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sig = ifac*sig;
	else
		sig = nfac*sig;
	end

	if hist
		dhist(iters) = dopt;
	end

	% Test convergence

	if sig < tol || dopt < tol
		converged = true;
		break
	end

end
