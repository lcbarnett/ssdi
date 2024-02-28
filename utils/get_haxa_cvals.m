function cvals = get_haxa_cvals(n,mdim,slev,datadir,verb)

% Get critical values for distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% mdim       vector of hyperplane dimensions
% slev       vector of significance levels
% datadir    stats file directory
% verb       verbose?
%
% cvals      critical angles (radians, in [0,pi/2] - i.e., unnormalised)

global local_data_dir

% Parameter defaults

if nargin < 2 || isempty(mdim),    mdim    = 1:n-1;          end % all
if nargin < 3 || isempty(slev),    slev    = [0.05 0.95];    end % standard tails
if nargin < 4 || isempty(datadir), datadir = local_data_dir; end % local data directory
if nargin < 5 || isempty(verb),    verb    = false;          end % not verbose

assert(isvector(mdim),'Hyperplane dimension(s) must be a vector');
assert(isvector(slev),'Significance level(s) must be a vector');

% Read stats file

assert(~isempty(datadir),'Stats directory path empty!');
if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
haxa_stats_file = fullfile(datadir,'haxa_stats.mat');
assert(exist(haxa_stats_file) == 2,'Hyperplane angle stats file ''%s'' not found',haxa_stats_file);
if verb, fprintf('Reading hyperplane angle stats from ''%s''... ',haxa_stats_file); end
load(haxa_stats_file);
if verb, fprintf('max. dimension = %d, sample size = %d\n',nmax,N); end

assert(n >= 2 && n <= nmax,        'Enclosing space dimension out of range 2 - %d',nmax);
assert(all(mdim > 0 & mdim < n),   'Some hyperplane dimension(s) out of range 1 - %d',n-1);
assert(all(slev >= 0 & slev <= 1), 'Some hyperplane significance levels(s) out of range 0 - 1');

haxa_cvaln = haxa_cval{n};

% Get critical values (interpolate critical values in haxa_cvaln for significance levels in haxa_slev)

nmdim = length(mdim);
nslev = length(slev);
cvals = nan(nmdim,nslev);
for k = 1:nslev
	for m = 1:nmdim
		cvals(m,k) = interp1(haxa_slev,haxa_cvaln(mdim(m),:),slev(k),'linear');
	end
end
