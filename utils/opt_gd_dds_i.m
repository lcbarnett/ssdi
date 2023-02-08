function [dopt,Qopt,converged,sig,iters,dhist] = opt_gd_dds_i(H,Qopt,maxiters,sig,gdls,tol,hist)

% Assumptions
%
% 1 - Qopt is orthonormal
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
	gtol = tol;
else
	stol = tol(1);
	dtol = tol(2);
	gtol = tol(3);
end

[n,m] = size(Qopt);

% Calculate dynamical dependence of initial projection

dopt  = trfun2dd(Q2L(Qopt),H);

if hist
	[~,g] = trfun2ddgrad(Qopt,H);  % dynamical dependence gradient
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dopt sig g];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% Move (hopefully) down gradient and orthonormalise

	[G,g] = trfun2ddgrad(Qopt,H);         % dynamical dependence gradient and magnitude
	Qtry  = orthonormalise(Qopt-sig*G/g); % gradient descent

	% Calculate dynamical dependence of trial projection

	dtry = trfun2dd(Q2L(Qtry),H);

	% If dynamical dependence smaller, accept move

	if dtry < dopt
		Qopt = Qtry;
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
