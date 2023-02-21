function [dd,L,converged,sig,iters,dhist] = opt_gd_ddx(CAK,L0,maxiters,sig,gdls,tol,hist)

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
	gtol = tol/10;
else
	stol = tol(1);
	dtol = tol(2);
	gtol = tol(3);
end

% Calculate proxy dynamical dependence of initial projection

L     = L0;
[G,g] = cak2ddxgrad(L,CAK); % proxydynamical dependence gradient and magnitude
dd    = cak2ddx(L,CAK);

if hist
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dd sig g];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters

	% Move (hopefully) down gradient and orthonormalise

	Ltry  = orthonormalise(L-sig*(G/g)); % gradient descent
	ddtry = cak2ddx(Ltry,CAK);

	% If dynamical dependence smaller, accept move and increase step size;
	% else reject move and decrease step size (similar to 1+1 ES)

	if ddtry < dd
		L     = Ltry;
		[G,g] = cak2ddxgrad(L,CAK); % proxy dynamical dependence gradient and magnitude
		dd    = ddtry;
		sig   = ifac*sig;
	else
		sig   = nfac*sig;
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
