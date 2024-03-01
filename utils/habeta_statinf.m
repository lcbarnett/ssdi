function [cval,pval,sig] = habeta_statinf(beta,n,m,slevel,tails,mhtc)

% Get critical values, p-values and significance for left- and
% right-tailed tests of the hyperplane axis angle beta statistic
% against the null hypothesis that the corresponding hyperplane
% is uniform random on the Grassmannian G(n,m). Under the null,
% the statistic is distributed as B(m/2,(n-m)/2); see habeta.m.
%
% A Beta statistic significantly close to 1 (right-tail) indicates
% significant participation of the channel corresponding to that
% axis with the hyperplane; a beta statistic significantly close
% to 0 (left-tail) indicates significant NON-participation.
%
% The significance level may optionally (should!) be adjusted for
% multiple hypotheses, using a Bonferroni correction.
%
% beta       Beta statistics (vector) = cos^2 of angles of hyperplane
%            with coordinate axes;  see metrics/habeta.m.
% n          dimension of enclosing Euclidean space
% m          dimension of hyperplane
% slevel     significance level (alpha)
% tails      'left' (non-participation), 'right' (participation) or 'both'
% mhtc       use Bonferroni multiple hypotheses test correction?
%
% If tails is set to 'both', pval and sig are returned as nbeta x 2
% matrices, where nbeta is the number of beta statistics supplied;
% the 1st column are the left-tail (non-participation) values, the
% 2nd column are the right-tail (participation) values; otherwise
% they are (column) vectors. cval is a (row) 2-vector (left- and
% right-tails).

if nargin < 4 || isempty(slevel), slevel = 0.05;   end
if nargin < 5 || isempty(tails),  tails  = 'both'; end
if nargin < 6 || isempty(mhtc),   mhtc   = true;   end

assert(isvector(beta),'Beta statistics must be a vector');
beta  = beta(:);      % ensure column vector
nbeta = length(beta); % number of stats to test

if mhtc
	slevel = slevel/(2*nbeta); % Bonferroni correction; 2*nbeta because left- and right-tailed tests!
end

% Beta distribution parameters

a = m/2;
b = (n-m)/2;
bcdf = betacdf(beta,a,b);

if     strcmpi(tails,'left')

	cval = betainv(slevel,a,b);
	pval = [bcdf 1-bcdf];

elseif strcmpi(tails,'right')

	cval = betainv(1-slevel,a,b);
	pval = [bcdf 1-bcdf];

elseif strcmpi(tails,'both')

	cval = betainv([slevel 1-slevel],a,b);
	pval = [bcdf 1-bcdf];

else
	error('Bad ''tails'' parameter');
end

sig  = pval <= slevel;
