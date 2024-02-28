function theta = make_haxa_stats(nmax,N,datadir)

% Calculate statistics (means, standard deviations, critical values) for
% distribution of angles of random hyperplanes with a coordinate axis.
%
% nmax       maximum enclosing space dimension
% N          stats sample size
% datadir    stats file directory
%
% Statistics are calculated for hyperplanes of dimension 1 .. n-1 from
% saved empirical distributions. Critical values are calculated for a set
% of 100 significance levels 0.0001 - 0.9999; linear interpolation may be
% used to estimate critical values for a given supplied significance level;
% see utils/get_haxa_cvals.m.

% Parameter defaults

if nargin < 2 || isempty(N),       N       = 100000;  end
if nargin < 3 || isempty(datadir), datadir = pwd;     end

% Significance levels

seg1 = 0.0001:0.0001:0.0004;
seg2 = 0.0005:0.0005:0.010;
seg3 = 0.015:0.005:0.10;
seg4 = 0.125:0.025:0.275;

haxa_slev = [0 seg1 seg2 seg3 seg4 1-fliplr(seg4) 1-fliplr(seg3) 1-fliplr(seg2) 1-fliplr(seg1) 1]; % 100 significance levels

assert(~isempty(datadir),'Stats directory path empty!');
if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator

haxa_mean = cell(nmax,1);
haxa_sdev = cell(nmax,1);
haxa_cval = cell(nmax,1);

for n = 2:nmax

	% Load hyperplane angle empirical distributions

	haxa_dist_file = fullfile(datadir,sprintf('haxa_dist_n%03d_N%d.mat',n,N));
	assert(exist(haxa_dist_file) == 2,'Hyperplane angle distribution file ''%s'' not found',haxa_dist_file);
	fprintf('Loading hyperplane angle distribution file: ''%s'' ... ',haxa_dist_file);
	load(haxa_dist_file);

	fprintf('calculating stats ... ');

	% Saved distributions are for m = 1 .. h = floor(n/2); reflect second half of distributions

	h = floor(n/2);
	theta(:,n-(1:h)) = pi/2-theta(:,1:h);

	% Calculate sample statistics

	haxa_mean{n} = mean(theta)'; % mean
	haxa_sdev{n} = std(theta)';  % std. deviation
	haxa_cval{n} = quantile(theta,haxa_slev)'; % critical values at significance levels corresponding to haxa_slev

	fprintf('done\n');
end

% Save results

haxa_stats_file = fullfile(datadir,'haxa_stats.mat');
fprintf('Saving hyperplane angle stats file: ''%s'' ... ',haxa_stats_file);
save(haxa_stats_file,'nmax','N','haxa_slev','haxa_mean','haxa_sdev','haxa_cval');
fprintf('done\n\n');
