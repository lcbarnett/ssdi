function [d,dd,converged,sig,i,L] = opt_ss1_es1p1_x(A,C,K,P0,iters,sig,ifac,nfac,tol)

[n,m] = size(P0);

% Calculate CAK sequence

CAK = CAK_seq(A,C,K);

% Orthonormalise initial projection

[L,M] = orthonormalise(P0);

% Calculate "fake" dynamical dependence of initial projection

dd = nan(iters,1);
d = ssdd(L,M,CAK);
dd(1) = d;

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
	dd(i) = d;

	% Test convergence

	if sig < tol || d < tol
		converged = true;
		break
	end

end
