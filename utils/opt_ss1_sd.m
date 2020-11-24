function [d,L] = opt_ss1_sd(A,C,K,P0,iters,sig)

[n,m] = size(P0);

% Calculate CAK sequence

CAK = CAK_seq(A,C,K);

% Orthonormalise initial projection

[L,M] = orthonormalise(P0);

% Calculate "fake" dynamical dependence of initial projection

d = ssdd1(L,M,CAK);

% Optimise

for i = 2:iters

	% "Mutate" projection and orthonormalise

	[Ltry,Mtry] = orthonormalise(L + sig*randn(n,m));

	% Calculate "fake" dynamical dependence of mutated projection

	dtry = ssdd1(Ltry,Mtry,CAK);

	% If dynamical dependence smaller, accept mutant

	if dtry < d
		L = Ltry;
		d = dtry;
	end

end
