function theta = subspaceb(A,B)

% An efficient subspace angle calulation for A, B orthonormal

if size(A,2) < size(B,2)
	P = B'*A;
else
	P = A'*B;
end
theta = acos(sqrt(min(eig(P'*P))));
