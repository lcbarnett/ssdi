function [LC,LK] = iss_perfect_dd(C,K,uniq)

% Generic perfect solutions exist only when 1 <= r < n.
% There is "overlap" between C and K solutions iff 2*r <= n.
%
% C solutions are unique iff m == n-r, K solutions iff m == r;
% otherwise an (arbitrary) selection of hierarchically-nested
% solutions is returned. If 'uniq' is set, then only the unique
% C and K solutions are returned.
%
% DDs should ALL be identically zero.

if nargin < 3 || isempty(uniq), uniq = true; end

[n,r] = size(C);
assert(r < n,'No generic perfect solutions if n >= r');

if uniq

	[LC,~] = orthonormalise([-C(r+1:n,:)/C(1:r,:) eye(n-r)]'); % m == n-r
	[~,LK] = orthonormalise([-K(:,1:r)\K(:,r+1:n);eye(n-r)] ); % m == r

else
	% C solutions (m = 1:n-r; solution is unique if m == n-r)

	LC = cell(n,1);
	for m = 1:n-r
		k = n-m;
		[LC{m},~] = orthonormalise([-C(k+1:n,:)/C(1:r,:) zeros(m,k-r) eye(m)]'); % the zero block fixes an arbitrary nested set of solutions

	end

	% K solutions (m = r:n-1; solution is unique if m == r)

	LK = cell(n,1);
	for m = r:n-1
		k = n-m;
		[~,LK{m}] = orthonormalise([-K(:,1:r)\K(:,m+1:n);zeros(m-r,k);eye(k)] ); % the zero block fixes an arbitrary nested set of solutions
	end
end
