function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_dd(A,C,K,Lopt,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

[n,m] = size(Lopt);

% Calculate dynamical dependence of initial projection

dopt = iss2dd(Lopt,A,C,K);

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

	% Calculate dynamical dependence of mutated projection

	dtry = iss2dd(Ltry,A,C,K);

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
