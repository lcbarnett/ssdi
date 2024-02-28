function cvals = get_haxa_cvals(n,mdim,slev,datadir)

% Get critical values for distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% mdim       vector of hyperplane dimensions
% slev       vector of significance levels
% datadir    stats file directory
%
% cvals      critical angles (radians, in [0,pi/2] - i.e., unnormalised)

global haxa_data_dir

% Parameter defaults

if nargin < 2 || isempty(mdim),    mdim    = 1:n-1;         end % all
if nargin < 3 || isempty(slev),    slev    = [0.05 0.95];   end % standard tails
if nargin < 4 || isempty(datadir), datadir = haxa_data_dir; end % local haxa stats directory

assert(isvector(mdim),'Hyperplane dimension(s) must be a vector');
assert(isvector(slev),'Significance level(s) must be a vector');

% Read stats file

assert(~isempty(datadir),'Stats directory path empty!');
if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
haxa_stats_file = fullfile(datadir,sprintf('haxa_stats_n%03d.mat',n));
assert(exist(haxa_stats_file) == 2,'Hyperplane angle stats file ''%s'' not found',haxa_stats_file);
fprintf('Reading hyperplane angle stats from ''%s''... ',haxa_stats_file);
load(haxa_stats_file);
fprintf('sample size = %d\n',N);

assert(all(slev > min(haxa_slev)-eps & slev < max(haxa_slev)-eps),'Some significance level(s) out of bounds (must be in range 0.001 - 0.999)');

nmdim = length(mdim);
nslev = length(slev);

% Get critical values (interpolate critical values in haxa_cval for significance levels in haxa_slev)

cvals = nan(nmdim,nslev);
for k = 1:nslev
	for m = 1:nmdim
		cvals(m,k) = interp1(haxa_slev,haxa_cval(mdim(m),:),slev(k),'linear');
	end
end
