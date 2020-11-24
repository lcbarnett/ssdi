function [d,L] = opt_ss_sd(A,C,K,P0,iters,sig)

[n,m] = size(P0);

% Orthonormalise initial projection

L = orthonormalise(P0);

% Calculate dynamical dependence of initial projection

d = ssdd(L,A,C,K);

% Optimise

for i = 2:iters

	% "Mutate" projection and orthonormalise

	Ltry = orthonormalise(L + sig*randn(n,m));

	% Calculate dynamical dependence of mutated projection

	dtry = ssdd(Ltry,A,C,K);

	% If dynamical dependence smaller, accept mutant

	if dtry < d
		L = Ltry;
		d = dtry;
	end

end
