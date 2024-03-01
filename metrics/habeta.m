function b = habeta(L)

% Returns squared cosines of angle between L (assumed orthonormal)
% and each of the coordinate axes ('Beta statistic').
%
% Under the null hypothesis that the hyperplane specified by the
% orthonormal n x m basis L is uniform random on the Grassmannian
% G(n,m), the per-channel statistics are marginally distributed as
% B(m/2,(n-m)/2); but note that they are not independent.

b = sum(L.^2,2);
