function [d,dd,converged,sig,i,L] = opt_ss_es1p1_x(A,C,K,P0,iters,sig,ifac,nfac,tol)

[n,m] = size(P0);

% Orthonormalise initial projection

L = orthonormalise(P0);

% Calculate dynamical dependence of initial projection

dd = nan(iters,1);
d = ssdd(L,A,C,K);
dd(1) = d;

% Optimise

converged = false;
for i = 2:iters

	% "Mutate" projection and orthonormalise

	Ltry = orthonormalise(L + sig*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = ssdd(Ltry,A,C,K);

	% If dynamical dependence smaller, accept mutant

	if dtry < d
		L = Ltry;
		d = dtry;
		sig = ifac*sig;
	else
		sig = nfac*sig;
	end
	dd(i) = d;

	% Test convergence

	if sig < tol || d < tol
		converged = true;
		break
	end

end
