function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_dds(H,Lopt,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

if isscalar(tol)
	stol = tol;
	dtol = tol;
else
	stol = tol(1);
	dtol = tol(2);
end

[n,m] = size(Lopt);

% Calculate dynamical dependence of initial projection

dopt = trfun2dd(Lopt,H);

if hist
	dhist = zeros(maxiters,2);
	dhist(1,:) = [dopt sig];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% "Mutate" projection and orthonormalise

	Ltry = orthonormalise(Lopt + sig*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = trfun2dd(Ltry,H);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sig  = ifac*sig;
	else
		sig  = nfac*sig;
	end

	if hist
		dhist(iters,:) = [dopt sig];
	end

	% Test convergence

	if sig < stol
		converged = 1;
		break
	end

	if dopt < dtol
		converged = 2;
		break
	end

end

if hist
	dhist = dhist(1:iters,:);
end
