function [F,err] = vou_to_dd(A,V,L,pstab)

% Calculate time-domain dynamical dependence rate of a specified
% linear subspace for a vector Ornstein-Uhlenbeck (VOU) process.
% Uses a state-space method which involves solving an associated
% continuous-time algebraic Riccati equation (CARE).
%
% A     - VOU coefficients matrix
% V     - VOU Wiener process covariance matrix
% L     - linear subspace basis
% pstab - stabilising perturbation for CARE (pstab = sqrt(eps) seems to work well)
%
% F     - DD rate for linear subspace
% err   - CARE error report number (zero if no error); run carerep(err) for eror message
%
% REFERENCES:
%
% (1) L. Barnett and A. K. Seth (2015): Granger causality for state-space models, Phys. Rev. E 91(4) Rapid Communication.
% (2) L. Barnett and A. K. Seth (2016): Detectability of Granger causality for subsampled continuous-time neurophysiological processes, J. Neurosci. Methods 275.
% (3) L. Barnett (2017): Granger causality rate for a vector Ornstein-Uhlenbeck process (working notes).
%
% (C) Lionel Barnett, May 2017

[n, n1]  = size(A); assert(n1 == n, 'VOU coeffcicients matrix must be square');
[n1,n2]  = size(V); assert(n1 == n2,'VOU covariance matrix must be square');
                    assert(n1 == n, 'VOU covariance matrix must be same size as coefficients matrix');

m = size(L,2);
if nargin > 3 && ~isempty(pstab)
	L = L + pstab*(rand(size(L))-0.5)/n;
end

[L,M] = orthonormalise(L);

A12 = L'*A*M;
A22 = M'*A*M;
V11 = L'*V*L;
V21 = M'*V*L;
V22 = M'*V*M;

[P22,~,~,rep] = icare(A22',A12',V22,V11,V21);
err = rep.Report;
if err ~= 0
	F = NaN;
	return
end

F = trace(V11\(A12*P22*A12'));
