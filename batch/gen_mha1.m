function gen_haxa_stats(n,N,C,datadir,dryrun)

% Calculate statistics (means, standard deviations, critical values) for
% distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% N          sample size
% C          number of "chunks" to divide sample into to avoid memory problems
% datadir    directory to save stats file
% dryrun     just calculate estimate of maximum chunk size and return
%
% Stats are calculated for hyperplanes of dimension 1 .. n-1. Critical values
% are calculated for a set of 50 significance levels as in the cointegration
% routine 'jcitest' (Econometrics Toolbox); interp1 (linear) may be used to
% interpolate critical values for a given supplied significance level alpha.

% Parameter defaults

if nargin < 3 || isempty(C),       C       = 10;      end
if nargin < 4 || isempty(datadir), datadir = tempdir; end
if nargin < 5 || isempty(dryrun),  dryrun  = false;   end

S = N/C; % C is number of chunks, S is chunksize
assert(C*S == N,'C must divide N exactly');

maxcMB = (S*(n-1)*n*8)/1000/1000;
fprintf('\nMaximum chunk size ~ %g MB\n',maxcMB);

if dryrun, return; end

% Sample angles for subspace dimension m = 1 .. n-1

theta = zeros(N,n-1);
m = 1;
fprintf('\nm = %2d of %2d\n',m,n-1);
for c = 1:C
	fprintf('\tchunk %2d of %2d\n',c,C);
	L = rand_orthonormal(n,m,S);
	theta((c-1)*S+1:c*S,m) = acos(abs(squeeze(L(1,:,:))));
end
clear L
for m = 2:n-1
	fprintf('\nm = %2d of %2d\n',m,n-1);
	for c = 1:C
		fprintf('\tchunk %2d of %2d\n',c,C);
		L = rand_orthonormal(n,m,S);
		theta((c-1)*S+1:c*S,m) = acos(sqrt(sum(squeeze(L(1,:,:)).^2)));
	end
	clear L
end

% Calculate sample statistics

q = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999]; % 50 quantiles as in 'jcitest'

haxa_mean = mean(theta)';       % mean
haxa_sdev = std(theta)';        % std. deviation
haxa_cval = quantile(theta,q)'; % critical values at significance levels corresponding to q

% Save results

if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
ffname = fullfile(datadir,sprintf('haxa_stats_n%03d_N%d.mat',n,N));
fprintf('\nSaving data file: ''%s'' ... ',ffname);
save(ffname,'n','N','haxa_mean','haxa_sdev','haxa_cval');
fprintf('done\n\n');
