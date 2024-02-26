function gen_haxa_dist(n,N,C,datadir,dryrun)

% Empirical distribution of angles of random hyperplanes with a coordinate axis.
%
% n          dimension of enclosing space
% N          stats sample size
% C          number of "chunks" to divide sample into to avoid memory problems
% datadir    distribution file directory
% dryrun     just calculate estimate of maximum chunk size and return
%
% Distributions are generated for hyperplanes of dimension 1 .. h = floor(n/2);
% dimensions up to n-1 may be obtained by relection (see make_haxa_stats.m). If
% dryrun = true, approximate memory footprint is calculated, and function exits.

% Parameter defaults

if nargin < 2 || isempty(N),       N       = 100000;  end
if nargin < 3 || isempty(C),       C       = 10;      end
if nargin < 4 || isempty(datadir), datadir = tempdir; end
if nargin < 5 || isempty(dryrun),  dryrun  = false;   end

S = N/C; % C is number of chunks, S is chunksize
assert(C*S == N,'C must divide N exactly');

n1 = n-1;
h = floor(n/2);

memreq = h*(N + n*S)*8;
if     memreq < 1000,           fprintf('\nMemory footprint ~ %d B\n',   memreq);
elseif memreq < 1000*1000,	    fprintf('\nMemory footprint ~ %.2f KB\n',memreq/1000);
elseif memreq < 1000*1000*1000, fprintf('\nMemory footprint ~ %.2f MB\n',memreq/1000/1000);
else,                           fprintf('\nMemory footprint ~ %.2f GB\n',memreq/1000/1000/1000);
end

if dryrun, return; end

% Sample angles for subspace dimension m = 1 .. h (only need first half - huge saving in computation!)

theta = zeros(N,h);
L = zeros(n,h,S); % preallocate
for m = h:-1:2
	fprintf('\nm = %2d of %2d\n',m,h);
	for c = 1:C
		fprintf('\tchunk %2d of %2d\n',c,C);
		L = rand_orthonormal(n,m,S);
		theta((c-1)*S+1:c*S,m) = acos(sqrt(sum(squeeze(L(1,:,:)).^2)));
	end
end
m = 1;
fprintf('\nm = %2d of %2d\n',m,h);
for c = 1:C
	fprintf('\tchunk %2d of %2d\n',c,C);
	L = rand_orthonormal(n,m,S);
	theta((c-1)*S+1:c*S,m) = acos(abs(squeeze(L(1,:,:))));
end

% Save results

if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
haxa_dist_file = fullfile(datadir,sprintf('haxa_dist_n%03d_N%d%s.mat',n,N));
fprintf('\nSaving hyperplane angle distribution file: ''%s'' ... ',haxa_dist_file);
save(haxa_dist_file,'n','N','theta');
fprintf('done\n\n');
