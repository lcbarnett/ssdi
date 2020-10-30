function HSV = make_nHSV(w,h)

if nargin < 2 || isempty(h), h = 240/360; end

n = size(w,1);

HSV = [h*ones(n,1) w ones(n,1)];
