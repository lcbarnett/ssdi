function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_dds(H,P0,maxiters,sig,ifac,nfac,tol,hist)

fres = size(H,3)-1;

[n,m] = size(P0);

% Orthonormalise initial projection

Lopt = orthonormalise(P0);

% Calculate dynamical dependence of initial projection

dopt = trfun2dd(Lopt,H);

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

	Ltry = orthonormalise(Lopt + sig*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = trfun2dd(Ltry,H);

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
