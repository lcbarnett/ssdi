function [dopt,Lopt,converged,sig,iters,dhist] = opt_gd_ddx(CAK,Lopt,maxiters,sig,gdls,tol,hist)

% Assumptions
%
% 1 - Lopt is orthonormal
% 2 - Residuals covariance matrix is identity

if isscalar(gdls)
	ifac = gdls;
	nfac = 1/ifac;
else
	ifac = gdls(1);
	nfac = gdls(2);
end

if isscalar(tol)
	stol = tol;
	dtol = tol;
	gtol = tol/10;
else
	stol = tol(1);
	dtol = tol(2);
	gtol = tol(3);
end

[n,m] = size(Lopt);

% Calculate proxy dynamical dependence of initial projection

dopt = cak2ddx(Lopt,CAK);

if hist
	[~,g] = cak2ddxgrad(Lopt,CAK);  % dynamical dependence gradient
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dopt sig g];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% Move (hopefully) down gradient and orthonormalise

	[G,g] = cak2ddxgrad(Lopt,CAK);        % dynamical dependence gradient and magnitude
	Ltry  = orthonormalise(Lopt-sig*G/g); % gradient descent

	% Calculate proxy dynamical dependence of mutated projection

	dtry = cak2ddx(Ltry,CAK);

	% If dynamical dependence smaller, accept mutant

	if dtry < dopt
		Lopt = Ltry;
		dopt = dtry;
		sig  = ifac*sig;
	else
		sig  = nfac*sig;
	end

	if hist
		dhist(iters,:) = [dopt sig g];
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

	if g < gtol
		converged = 3;
		break
	end

end

if hist
	dhist = dhist(1:iters,:);
end
