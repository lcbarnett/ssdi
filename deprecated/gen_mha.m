%-------------------------------------------------------------------------------

defvar('n',        40     );
defvar('N',        100    );
defvar('mangle',   true   );

%-------------------------------------------------------------------------------

theta = zeros(N,n-1,n-1);
k=1;
parfor k=1:(n-1) % set parpool first (if you want)
	for m = 1:k
		fprintf('k = %2d of %2d, m = %2d\n',k,n-1,m);
		Lk = [eye(k);zeros(n-k,k)];
		L = rand_orthonormal(n,m,N);
		for i = 1:N
			theta(i,m,k) = gmetric(L(:,:,i),Lk,mangle);
		end
	end
end
thetam = symmetrise(squeeze(mean(theta))); % mean
thetas = symmetrise(squeeze(std(theta)));  % std.deviation
thetae = thetas/sqrt(N);                   % std. error

global PMDATAROOT
if mangle
	matfile = sprintf('mhax_N%d_n%03d.mat',N,n);
else
	matfile = sprintf('mham_N%d_n%3d.mat',N,n);
end
ffname = fullfile(PMDATAROOT,'metadata','mha',matfile);
fprintf('Saving data file: ''%s'' ... ',ffname);
save(ffname,'N','n','mangle','thetam','thetas','thetae');
fprintf('done\n\n');
