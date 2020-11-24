function [d,converged,sig,i,L] = opt_ss1_es1p1(A,C,K,P0,iters,sig,ifac,nfac,tol)

[n,m] = size(P0);

% Calculate CAK sequence

CAK = CAK_seq(A,C,K);

% Orthonormalise initial projection

[L,M] = orthonormalise(P0);

% Calculate "fake" dynamical dependence of initial projection

d = ssdd1(L,M,CAK);

% Optimise

converged = false;
for i = 2:iters

	% "Mutate" projection and orthonormalise

	[Ltry,Mtry] = orthonormalise(L + sig*randn(n,m));

	% Calculate "fake" dynamical dependence of mutated projection

	dtry = ssdd1(Ltry,Mtry,CAK);

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
