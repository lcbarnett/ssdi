function [dopt,Lopt,converged,sig,iters,dhist] = opt_es_sldd(A,C,K,Lopt,maxiters,sig,ifac,nfac,tol,hist)

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

% Calculate 'ldwork' parameter for slidare

[~,info] = slidare(A',C'*L,K*K',[],(K*L)');
ldwork = info.ldwork;

% Calculate dynamical dependence of initial projection

dopt = sliss2dd(Lopt,A,C,K,ldwork);

if hist
	dhist = zeros(maxiters,3);
	i = 1; dhist(i,:) = [1 dopt sig];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% "Mutate" projection and orthonormalise

	Ltry = orthonormalise(Lopt + sig*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = sliss2dd(Ltry,A,C,K,ldwork);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		if hist
			i = i+1; dhist(i,:) = [iters dopt sig];
			dopt = dtry;
			sig  = ifac*sig;
			i = i+1; dhist(i,:) = [iters dopt sig];
		else
			dopt = dtry;
			sig  = ifac*sig;
		end
	else
		sig = nfac*sig;
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
	i = i+1; dhist(i,:) = [iters dopt sig];
	dhist = dhist(1:i,:);
end
