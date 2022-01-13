function [dd,L,converged,sig,iters,dhist] = opt_gd_dds(H,L,maxiters,sig,ifac,nfac,tol,hist)

% Assumptions
%
% 1 - L is orthonormal
% 2 - Residuals covariance matrix is identity

if isscalar(tol)
	stol = tol;
	dtol = tol;
	gtol = tol;
else
	stol = tol(1);
	dtol = tol(2);
	gtol = tol(3);
end

[n,m] = size(L);

% Calculate dynamical dependence of initial projection

dd    = trfun2dd(L,H);         % current dynamical dependence
grad  = trfun2ddgrad(L,H);     % dynamical dependence gradient
mgrad = sqrt(sum(grad(:).^2)); % magnitude of gradient vector

if hist
	dhist = zeros(maxiters,3);
	dhist(1,:) = [dd sig mgrad];
else
	dhist = [];
end

% Optimise

converged = 0;
for iters = 2:maxiters


	% Move (hopefully) down gradient and orthonormalise

	grad  = trfun2ddgrad(L,H);            % dynamical dependence gradient
	mgrad = sqrt(sum(grad(:).^2));        % magnitude of gradient vector
	L = orthonormalise(L-sig*grad/mgrad); % gradient descent
	ddold = dd;                           % save previous dynamical dependence
	dd = trfun2dd(L,H);                   % dynamical dependence of new projection

	% Adjust learning rate
%{
	if dd < ddold
		sig  = ifac*sig;
	else
		sig  = nfac*sig;
	end
%}
	if hist
		dhist(iters,:) = [dd sig mgrad];
	end

	% Test convergence

	if sig < stol
		converged = 1;
		break
	end

	if dd < dtol
		converged = 2;
		break
	end

	if mgrad < gtol
		converged = 3;
		break
	end

end

if hist
	dhist = dhist(1:iters,:);
end
