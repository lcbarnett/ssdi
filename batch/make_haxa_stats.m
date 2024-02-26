function theta = make_haxa_stats(n,N,datadir,hstag)

% Calculate statistics (means, standard deviations, critical values) for
% distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% N          stats sample size
% datadir    stats file directory
% hstag      stats file ID tag
%
% Statistics are calculated for hyperplanes of dimension 1 .. n-1 from saved
% empirical distributions. Critical values are calculated for a set of 50
% significance levels as in the cointegration routine 'jcitest' (Econometrics
% Toolbox); interp1 (linear) may be used to interpolate critical values for
% a given supplied significance level.

% Parameter defaults

if nargin < 2 || isempty(N),       N       = 100000;  end
if nargin < 3 || isempty(datadir), datadir = tempdir; end
if nargin < 4                      hstag   = '';      end % optional stats ID tag

% Load hyperplane angle empirical distributions

if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
haxa_dist_file = fullfile(datadir,sprintf('haxa_dist_n%03d_N%d.mat',n,N));
assert(exist(haxa_dist_file) == 2,'Hyperplane angle distribution file ''%s'' not found',haxa_dist_file);
fprintf('\nLoading hyperplane angle distribution file: ''%s'' ... ',haxa_dist_file);
load(haxa_dist_file);
fprintf('done\n\n');

fprintf('Calculating stats ... ');

% Saved distributions are for m = 1 .. h = floor(n/2); reflect second half of distributions

h = floor(n/2);
theta(:,n-(1:h)) = pi/2-theta(:,1:h);

% Calculate sample statistics

haxa_slev = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999]; % 50 significance levels as in 'jcitest'
haxa_mean = mean(theta)'; % mean
haxa_sdev = std(theta)';  % std. deviation
haxa_cval = quantile(theta,haxa_slev)'; % critical values at significance levels corresponding to haxa_slev

fprintf('done\n\n');

% Save results

if ~isempty(hstag), hstag = ['_' hstag]; end
haxa_stats_file = fullfile(datadir,sprintf('haxa_stats_n%03d%s.mat',n,hstag));
fprintf('Saving hyperplane angle stats file: ''%s'' ... ',haxa_stats_file);
save(haxa_stats_file,'n','N','haxa_mean','haxa_sdev','haxa_cval','haxa_slev');
fprintf('done\n\n');
