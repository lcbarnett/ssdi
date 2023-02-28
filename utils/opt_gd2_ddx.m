function [dd,L,converged,sig,iters,dhist] = opt_gd2_ddx(CAK,L0,maxiters,gdsig0,gdls,tol,hist)

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
sig   = gdsig0;

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

	L     = orthonormalise(L-sig*(G/g)); % gradient descent
	[G,g] = cak2ddxgrad(L,CAK); % proxy dynamical dependence gradient and magnitude
	ddnew = cak2ddx(L,CAK);

	% If dynamical dependence smaller, accept move and increase step size;
	% else reject move and decrease step size (similar to 1+1 ES)

	if ddnew < dd
		dd  = ddnew;
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
