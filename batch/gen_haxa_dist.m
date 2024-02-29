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
if nargin < 4 || isempty(datadir), datadir = pwd;     end
if nargin < 5 || isempty(dryrun),  dryrun  = false;   end

S = N/C; % C is number of chunks, S is chunksize
assert(C*S == N,'C must divide N exactly');

h = floor(n/2);

memreq = (h*N+n*S)*8;
if     memreq < 1000,           fprintf('\nMemory footprint ~ %d B\n',   memreq);
elseif memreq < 1000*1000,	    fprintf('\nMemory footprint ~ %.2f KB\n',memreq/1000);
elseif memreq < 1000*1000*1000, fprintf('\nMemory footprint ~ %.2f MB\n',memreq/1000/1000);
else,                           fprintf('\nMemory footprint ~ %.2f GB\n',memreq/1000/1000/1000);
end

if dryrun, return; end

% Sample angles for subspace dimension m = 1 .. h (only need first half)

theta = zeros(N,h); % pre-allocate
v     = zeros(n,S); % pre-allocate

m = 1;
printf('m = %2d of %2d\n',m,h);
for c = 1:C
	fprintf('\tchunk %2d of %2d\n',c,C);
	v = randn(n,S);
	theta((c-1)*S+1:c*S,m) = acos(sqrt((v(1,:).^2)./sum(v.^2))); % angles with 1st axis
end
for m = 2:h
	fprintf('m = %2d of %2d\n',m,h);
	for c = 1:C
		fprintf('\tchunk %2d of %2d\n',c,C);
		v = randn(n,S);
		theta((c-1)*S+1:c*S,m) = acos(sqrt(sum(v(1:m,:).^2)./sum(v.^2))); % angles with hyperplane formed by 1st m axes
	end
end

% Save results

if datadir(end) == filesep, datadir = datadir(1:end-1); end % strip trailing file path separator
haxa_dist_file = fullfile(datadir,sprintf('haxa_dist_n%03d_N%d%s.mat',n,N));
fprintf('\nSaving hyperplane angle distribution file: ''%s'' ... ',haxa_dist_file);
save(haxa_dist_file,'n','N','theta');
fprintf('done\n\n');
