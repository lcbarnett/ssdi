function cvals = get_haxa_cvals(n,mdim,slev,haxa_stats_file)

% Get critical values for distribution of angles of random hyperplanes with a coordinate axis.

% Read stats file

assert(exist(haxa_stats_file,'file'),'Angle stats file ''%s'' not found',haxa_stats_file);
load(haxa_stats_file);
fprintf('Read angle stats from ''%s''\n',haxa_stats_file);

assert(isvector(mdim),'Hyperplane dimension(s) must be a vector');
assert(isvector(slev),'Significance level(s) must be a vector');
assert(all(slev > min(haxa_slev)-eps/2 & slev < max(haxa_slev)-eps/2),'Some significance level(s) out of bounds');

nmdim = length(mdim);
nslev = length(slev);

% Get critical values (interpolate critical values in haxa_cval for significance levels in haxa_slev)

cvals = nan(nmdim,nslev);
for k = 1:nslev
	for m = 1:nmdim
		cvals(m,k) = interp1(haxa_slev,haxa_cval(mdim(m),:),slev(k),'linear');
	end
end
