function theta = gen_haxa_stats(n,N,C,datadir,hstag,dryrun)

% Calculate statistics (means, standard deviations, critical values) for
% distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% N          stats sample size
% C          number of "chunks" to divide sample into to avoid memory problems
% datadir    stats file directory
% hstag      stats file ID tag
% dryrun     just calculate estimate of maximum chunk size and return
%
% Stats are calculated for hyperplanes of dimension 1 .. n-1. Critical values
% are calculated for a set of 50 significance levels as in the cointegration
% routine 'jcitest' (Econometrics Toolbox); interp1 (linear) may be used to
% interpolate critical values for a given supplied significance level alpha.

% Parameter defaults

if nargin < 2 || isempty(N),       N       = 100000;  end
if nargin < 3 || isempty(C),       C       = 10;      end
if nargin < 4 || isempty(datadir), datadir = tempdir; end
if nargin < 5                      hstag   = '';      end % no tag
if nargin < 6 || isempty(dryrun),  dryrun  = false;   end

S = N/C; % C is number of chunks, S is chunksize
assert(C*S == N,'C must divide N exactly');

n1 = n-1;

memreq = n1*(N + n*S + 2 + 50)*8/1000/1000;
fprintf('\nMemory requirements: at least %g MB\n',memreq);

if dryrun, return; end

% Sample angles for subspace dimension m = 1 .. n-1

theta = zeros(N,n1);
L = zeros(n,n1,S); % preallocate
for m = n1:-1:2
	fprintf('\nm = %2d of %2d\n',m,n1);
	for c = 1:C
		fprintf('\tchunk %2d of %2d\n',c,C);
		L = rand_orthonormal(n,m,S);
		theta((c-1)*S+1:c*S,m) = acos(sqrt(sum(squeeze(L(1,:,:)).^2)));
	end
end
m = 1;
fprintf('\nm = %2d of %2d\n',m,n1);
for c = 1:C
	fprintf('\tchunk %2d of %2d\n',c,C);
	L = rand_orthonormal(n,m,S);
	theta((c-1)*S+1:c*S,m) = acos(abs(squeeze(L(1,:,:))));
end

% Calculate sample statistics

haxa_slev = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999]; % 50 significance levels as in 'jcitest'
haxa_mean = mean(theta)'; % mean
haxa_sdev = std(theta)';  % std. deviation
haxa_cval = quantile(theta,haxa_slev)'; % critical values at significance levels corresponding to haxa_slev

% Save results

if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
if ~isempty(hstag), hstag = ['_' hstag]; end
haxa_stats_file = fullfile(datadir,sprintf('haxa_stats_n%03d%s.mat',n,hstag));
fprintf('\nSaving data file: ''%s'' ... ',haxa_stats_file);
save(haxa_stats_file,'n','N','haxa_mean','haxa_sdev','haxa_cval','haxa_slev');
fprintf('done\n\n');
