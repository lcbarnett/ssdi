function [ab,B] = habetax(n,m1,m2,N,C,fignum,dryrun)

% Returns fitted Beta parameters for the distribution  of the
% squared cosine of the angle between uniform random hyperplanes
% of dimensions m1, m2 in an enclosing space of dimension n, based
% on an empirical sample of size S. Optionally returns the sample
% itself in B.
%
% NOTE: Unless m1 = 1, the Beta distribution is approximate.
%
% The computation is segmented into C "chunks" of size N/C to avoid
% memory allocation issues.

assert(m1 > 0 && m2 < n,'Hyperplane dimensions out of bounds');
assert(m1 <= m2,'First hyperplane dimension cannot be larger than second');

% Parameter defaults

if nargin < 4 || isempty(N),       N       = 100000;  end
if nargin < 5 || isempty(C),       C       = 10;      end
if nargin < 6 || isempty(fignum),  fignum  = 1;       end
if nargin < 7 || isempty(dryrun),  dryrun  = 0;       end

S = N/C; % C is number of chunks, S is chunksize
assert(C*S == N,'C must divide N exactly');

ab = [];
B = [];

memreq = (C+n*m1)*S*8;
if     memreq < 1000,           fprintf('\nMemory footprint ≈ %d B\n\n',   memreq);
elseif memreq < 1000*1000,	    fprintf('\nMemory footprint ≈ %.2f KB\n\n',memreq/1000);
elseif memreq < 1000*1000*1000, fprintf('\nMemory footprint ≈ %.2f MB\n\n',memreq/1000/1000);
else,                           fprintf('\nMemory footprint ≈ %.2f GB\n\n',memreq/1000/1000/1000);
end

if isunix % Windos, MacOS - you're welcome to submit a patch
	[~,ma] = system('free -b | awk ''NR==2 {print $7}''');
	memavail = str2num(ma);
	if memreq > (3/4)*memavail
		fprintf('That''s pushing it a bit: only ≈ %.2f GB available. Try specifying more "chunks" (''C'' parameter)\n\n',memavail/1000/1000/1000);
		if dryrun >= 0
			fprintf(2,'To force execution anyway, set ''dryrun'' parameter to something negative\n\n');
			return
		end
	end
end

if dryrun > 0, return; end

B = zeros(S,C);
X = zeros(n,m1,S);
for c = 1:C
	fprintf('chunk %3d of %3d\n',c,C);
	X = rand_orthonormal(n,m1,S);
	for s = 1:S
		P = X(1:m2,:,s);
		B(s,c) = min(eig(P'*P,'nobalance'));
	end
end
B = reshape(B,N,1);
ab = betafit(B);

if fignum > 0
	x = linspace(0,1,1001)';
	bpdf = betapdf(x,ab(1),ab(2));
	figure(fignum);
	clf;
	histogram(B,'Normalization','pdf');
	hold on
	plot(x,bpdf,'LineWidth',1.5);
	hold off
	xlim([0,1]);
	xlabel('Beta');
	ylabel('PDF');
	title(sprintf('n = %d, m1 = %d, m2 = %d\n\nB(%.2f, %.2f)\n',n,m1,m2,ab(1),ab(2)));
end
