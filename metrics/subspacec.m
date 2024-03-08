function theta = subspacec(A,B)

% An efficient subspace angle calulation for A, B orthonormal

if size(A,2) < size(B,2)
	theta = asin(norm(A - B*(B'*A)));
else
	theta = asin(norm(B - A*(A'*B)));
end
