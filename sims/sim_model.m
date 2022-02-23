%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up SS or VAR model, calculate transfer function, CAK sequence, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',  tempdir    ); % model directory
defvar('modname', 'sim_model'); % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varmod = exist('G','var'); % for a VAR model, specify a connectivity matrix or a scalar dimension
if varmod % VAR model
	if isscalar(G), n = G; G = ones(n); else, n = size(G,1); end;
	defvar('r', 7   ); % VAR model order
	defvar('w', 1   ); % VAR coefficients decay parameter
else      % fully-connected state-space model
	defvar('n', 9   ); % microscopic dimension
	defvar('r', 3*n ); % hidden state dimension
end

defvar('rho',     0.9   ); % spectral norm (< 1)
defvar('rmii',    1     ); % residuals multiinformation; 0 for zero correlation
defvar('fres',    []    ); % frequency resolution (empty for automatic)
defvar('nsics',   0     ); % number of samples for spectral integration check (0 for no check)
defvar('mseed',   0     ); % model random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gvprog',  'neato' ); % GraphViz program/format (also try 'neato', 'fdp')
defvar('gvdisp',  true    ); % GraphViz display? Empty for no action, true to display, false to just generate files)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random VAR or ISS model
rstate = rng_seed(mseed);
V0 = corr_rand(n,rmii); % residuals covariance matrix
if varmod
	ARA0 = var_rand(G,r,rho,w);             % random VAR model
	gc = var_to_pwcgc(ARA0,V0);             % causal graph
	[ARA,V] = transform_var(ARA0,V0);       % transform model to decorrelated-residuals form
	[A,C,K] = var_to_ss(ARA);               % equivalent ISS model
	if isempty(fres)
		[fres,ierr] = var2fres(ARA,V);
	end
	CAK = ARA;                              % CAK sequence for pre-optimisation
	H = var2trfun(ARA,fres);                % transfer function
	mdescript = sprintf('%d-variable VAR(%d)',n,r);
else
	[A0,C0,K0] = iss_rand(n,r,rho);         % random ISS model
	gc = ss_to_pwcgc(A0,C0,K0,V0);          % causal graph
	[A,C,K,V] = transform_ss(A0,C0,K0,V0);  % transform model to decorrelated-residuals form
	if isempty(fres)
		[fres,ierr] = ss2fres(A,C,K,V);
	end
	CAK = iss2cak(A,C,K);                   % CAK sequence for pre-optimisation
	H = ss2trfun(A,C,K,fres);               % transfer function
	mdescript = sprintf('%d-variable ISS(%d)',n,r);
end
fprintf('\nFrequency resolution = %d (integration error = %e)\n\n',fres,ierr);
rng_restore(rstate);

% Spectral resolution accuracy check

if nsics > 0
	assert(exist('mdim','var'),'For spectral accuracy check, must supply macro dimension ''mdim''');
	fprintf('Spectral DD accuracy check (frequency resolution = %d) ... ',fres);
	derr = dds_check(A,C,K,H,m,nsics); % spectral integration check
	fprintf('integration error = %e\n\n',derr);
	if derr > 1e-12, fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n'); end
end

% Model info

fprintf('--------------------------------------------\n');
fprintf('Model                : %s\n',mdescript);
fprintf('--------------------------------------------\n');
fprintf('Dimension            : %d\n',n);
fprintf('Complexity (CAK)     : %d x %d x %d\n',size(CAK,1),size(CAK,2),size(CAK,3));
fprintf('Frequency resolution : %d\n',size(H,3)-1);
fprintf('--------------------------------------------\n\n');

% Optionally display causal graph

if ~isempty(gvdisp)
	eweight = gc/nanmax(gc(:));
	gfile = fullfile(tempdir,'sim_model');
	wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);
	fprintf('\n');
end

% Save model

modfile = [fullfile(moddir,modname) '.mat'];
fprintf('*** saving model in ''%s''... ',modfile);
save(modfile,'V0','CAK','H','mdescript');
fprintf('done\n\n');
