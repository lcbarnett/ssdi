function HSV = make_cHSV(w,h)

if nargin < 2 || isempty(h), h = 0/360; end

n = size(w,1);

HSV = nan(n,n,3);
HSV(:,:,1) = h*ones(n);
HSV(:,:,2) = w;
HSV(:,:,3) = ones(n);
