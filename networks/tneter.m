function G = tneter(n,p)

% Erdos-Renyi random (directed) graph with self-connections

G = eye(n);
G(rand(n) < p) = 1;
