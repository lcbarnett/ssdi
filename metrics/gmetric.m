function d = gmetric(L1,L2,maxangle)

% Metric on the Grassmanian manifold, normalised to lie in [0,1]

if nargin < 3 || isempty(maxangle), maxangle = true; end

theta = subspacea(L1,L2);
if maxangle
	d = max(theta)/(pi/2);
else
	d = sqrt(mean(theta.^2))/(pi/2); % default
end
