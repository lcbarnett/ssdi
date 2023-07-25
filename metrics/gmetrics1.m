function d = gmetrics1(L)

% Returns normalised angle between L and each of the coordinate axes

d = acos(sqrt(sum(L.^2,2)))/(pi/2);
