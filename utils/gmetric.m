function d = gmetric(L1,L2,maxangle)

% Metric on the Grassmanian manifold, normalised to lie in [0,1]

theta = subspacea(L1,L2);
if nargin > 2 && maxangle
	d = max(theta)/(pi/2);
else
	d = sqrt(mean(theta.^2))/(pi/2); % default
end
