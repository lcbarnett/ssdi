function [dd,L,converged,sig,iters,dhist] = opt_gd1_dds(H,L0,maxiters,gdsig0,gdls,tol,hist)

% Assumptions
%
% 1 - L is orthonormal
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

% Calculate dynamical dependence of initial projection

L     = L0;
[G,g] = trfun2ddgrad(L,H); % dynamical dependence gradient and magnitude
dd    = trfun2dd(L,H);
sig   = gdsig0;

if hist
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dd sig g];
else
	dhist = [];
end

% Optimise: gradient descent

converged = 0;
for iters = 2:maxiters

	% Move (hopefully) down gradient and orthonormalise

	Ltry  = orthonormalise(L-sig*(G/g)); % gradient descent
	ddtry = trfun2dd(Ltry,H);

	% If dynamical dependence smaller, accept move and increase step size;
	% else reject move and decrease step size (similar to 1+1 ES)

	if ddtry < dd
		L     = Ltry;
		[G,g] = trfun2ddgrad(L,H); % dynamical dependence gradient and magnitude
		dd  = ddtry;
		sig = ifac*sig;
	else
		sig = nfac*sig;
	end

	if hist
		dhist(iters,:) = [dd sig g];
	end

	% Test convergence

	if     sig < stol
		converged = 1;
		break
	elseif dd < dtol
		converged = 2;
		break
	elseif g  < gtol
		converged = 3;
		break
	end

end

if hist
	dhist = dhist(1:iters,:);
end
