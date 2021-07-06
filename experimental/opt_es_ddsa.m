function [dopt,Lopt,converged,sopt,iters,dhist] = opt_es_ddsa(H,Lopt,maxiters,sopt,ssig,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

[n,m] = size(Lopt);

% Calculate dynamical dependence of initial projection

dopt = trfun2dd(Lopt,H);

if hist
	dhist = nan(maxiters,2);
	dhist(1) = dopt;
	dhist(2) = sopt;
else
	dhist = [];
end

% Optimise

hssig2 = (ssig*ssig)/2;

converged = false;
for iters = 2:maxiters

	% "Mutate" projection and orthonormalise

	stry = sopt*exp(ssig*randn - hssig2); % log-normal
	Ltry = orthonormalise(Lopt + sopt*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = trfun2dd(Ltry,H);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sopt = stry;
	end

	if hist
		dhist(iters,1) = dopt;
		dhist(iters,2) = sopt;
	end

	% Test convergence

	if sopt < tol || dopt < tol
		converged = true;
		break
	end

end
