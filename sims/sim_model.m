%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up SS or VAR model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varmod = exist('CON','var'); % for a VAR model, specify a connectivity matrix or a scalar dimension
if varmod % VAR model
	if isscalar(CON), n = CON; CON = ones(n); else, n = size(CON,1); end;
	defvar('r', 7   ); % VAR model order
	defvar('w', 1   ); % VAR coefficients decay parameter
else      % fully-connected state-space model
	defvar('n', 9   ); % microscopic dimension
	defvar('r', 3*n ); % hidden state dimension
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('moddir',   tempdir    ); % model directory
defvar('modname', 'sim_model' ); % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('rho',      0.9   ); % spectral norm (< 1)
defvar('rmii',     1     ); % residuals multiinformation; 0 for zero correlation
defvar('fres',     []    ); % frequency resolution (empty for automatic)
defvar('nsics',    0     ); % number of samples for spectral integration check (0 for no check)
defvar('mseed',    0     ); % model random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('gvprog',  'neato' ); % GraphViz program/format (also try 'neato', 'fdp')
defvar('gvdisp',   true   ); % GraphViz display? Empty for no action, true to display, false to just generate files)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Generate random VAR or ISS model
%
% NOTE: parameters, etc., with 0 suffix are for untransformed model!

rstate = rng_seed(mseed);
V0 = corr_rand(n,rmii); % residuals covariance (actually correlation) matrix
if varmod
	mdescript = sprintf('%d-variable VAR(%d)',n,r);
	ARA0 = var_rand(CON,r,rho,w);           % random VAR model
	gc = var_to_pwcgc(ARA0,V0);             % causal graph
	[A0,C0,K0] = var_to_ss(ARA0);           % equivalent ISS model
	if isempty(fres)
		[fres,ierr] = var2fres(ARA0,V0);
	end
	modcomp = r*n*n + (n*(n+1))/2;          % model complexity
else
	mdescript = sprintf('%d-variable ISS(%d)',n,r);
	[A0,C0,K0] = iss_rand(n,r,rho);         % random ISS model
	gc = ss_to_pwcgc(A0,C0,K0,V0);          % causal graph
	if isempty(fres)
		[fres,ierr] = ss2fres(A0,C0,K0,V0);
	end
	modcomp = (2*n+1)*r + (n*(n+1))/2;      % model complexity
end
rng_restore(rstate);

% Spectral resolution accuracy check

if nsics > 0
	assert(exist('mdim','var'),'For spectral accuracy check, must supply macro dimension ''mdim''');
	derr = dds_check(A,C,K,H,mdim,nsics); % spectral integration check
	if derr > 1e-12, fprintf(2,'WARNING: spectral DD calculation may be inaccurate!\n\n'); end
end

% Model info

fprintf('\n---------------------------------------\n');
fprintf('Model            : %s\n',mdescript);
fprintf('---------------------------------------\n');
fprintf('Dimension        : %d\n',n);
fprintf('Complexity       : %d\n',modcomp);
fprintf('Freq. resolution : %d\n',fres);
fprintf('---------------------------------------\n\n');

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
if varmod
	save(modfile,'mdescript','varmod','n','r','rho','rmii','V0','ARA0','A0','C0','K0','gc','fres','modcomp');
else
	save(modfile,'mdescript','varmod','n','r','rho','rmii','V0',       'A0','C0','K0','gc','fres','modcomp');
end
fprintf('done\n\n');
