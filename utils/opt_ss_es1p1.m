function [d,converged,sig,i,L] = opt_ss_es1p1(A,C,K,L,iters,sig,ifac,nfac,tol)

[n,m] = size(L);

% Calculate dynamical dependence of initial projection

d = ssdd(L,A,C,K);

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

	% Test convergence

	if sig < tol || d < tol
		converged = true;
		break
	end

end
