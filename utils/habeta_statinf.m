function [cval,pval,sig] = habeta_statinf(beta,n,m,slevel,mhtc)

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
% beta       beta statistics (vector) = cos^2 of angles of hyperplane
%            with coordinate axes;  see habeta.m.
% n          dimension of enclosing Euclidean space
% m          dimension of hyperplane
% slevel     significance level (alpha)
% mhtc       multiple-hypothesis test correction (see MVGC2 stats/significance.m)
%
% pval and sig are returned as nbeta x 2 matrices, where nbeta is
% the number of beta statistics supplied; the 1st column are the
% left-tail (non-participation) values, the 2nd column the right-
% tail (participation) values. cval is a 2-vector (for left- and
% right-tails).

if nargin < 4 || isempty(slevel), slevel = 0.05;   end
if nargin < 5 || isempty(mhtc),   mhtc   = 'NONE'; end

assert(isvector(beta),'Beta statistics must be a vector');
beta  = beta(:); % ensure column vector

a = m/2;
b = (n-m)/2;

cval = betainv([slevel 1-slevel],a,b); % left-tail, right-tail

if nargout > 1

	bcdf = betacdf(beta,a,b);
	pval = [bcdf 1-bcdf];  % left-tail, right-tail

	if nargout > 2
		sig = significance(pval,slevel,mhtc);
	end
end
