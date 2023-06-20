function [LC,LK] = iss_perfect_dd(C,K)

% Generic perfect solutions exist only when 1 <= r < n.
% There is "overlap" between C and K solutions iff 2*r <= n.
%
% DDs should ALL be zero.

[n,r] = size(C);
assert(r < n,'No generic perfect solutions if n >= r');

% C solutions (m = 1:n-r; solution is UNIQUE if m == n-r)

LC = cell(n,1);
for m = 1:n-r
	LC{m} = orthonormalise([-C(n-m+1:n,:)/C(1:n-m,:) eye(m)]');
end

% K solutions (m = r:n-1; solution is UNIQUE if m == r)

LK = cell(n,1);
for m = r:n-1
	[~,LK{m}] = orthonormalise([-K(:,1:m)\K(:,m+1:n);eye(n-m)]);
end
