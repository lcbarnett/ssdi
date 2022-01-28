%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up SS or VAR model, calculate transfer function, CAK sequence, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('moddir',   'var'), moddir   = tempdir;    end % model directory
if ~exist('fnamem',   'var'), fnamem   = 'model_dd'; end % model filename root

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varmod = exist('G','var'); % for a VAR model, specify a connectivity matrix or a scalar dimension
if varmod % VAR model
	if isscalar(G), n = G; G = ones(n); else, n = size(G,1); end;
	if ~exist('r',    'var'), r = 7;                 end % VAR model order
	if ~exist('w',    'var'), w = 1;                 end % VAR coefficients decay parameter
else      % fully-connected state-space model
	if ~exist('n',    'var'), n = 9;                 end % microscopic dimension
	if ~exist('r',    'var'), r = 3*n;               end % hidden state dimension
end

if ~exist('rho',      'var'), rho      = 0.9;        end % spectral norm (< 1)
if ~exist('rmii',     'var'), rmii     = 1;          end % residuals multiinformation; 0 for zero correlation
if ~exist('fres',     'var'), fres     = [];         end % frequency resolution (empty for automatic)
if ~exist('sitol',    'var'), sitol    = 1e-12;      end % spectral integration tolerance
if ~exist('siminp2',  'var'), siminp2  = 6;          end % spectral integration freq. res. min power of 2
if ~exist('simaxp2',  'var'), simaxp2  = 14;         end % spectral integration freq. res. max power of 2
if ~exist('nsics',    'var'), nsics    = 0;          end % number of samples for spectral integration check (0 for no check)
if ~exist('mseed',    'var'), mseed    = 0;          end % model random seed (0 to use current rng state)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('gvprog',   'var'), gvprog   = 'neato'; end % GraphViz program/format (also try 'neato', 'fdp')
if ~exist('gvdisp',   'var'), gvdisp   = true;    end % GraphViz display? Empty for no action, true to display, false to just generate files)

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
		[fres,ierr] = var2fres(ARA,V,[sitol,siminp2,simaxp2]);
		if isnan(fres) % failed!
			fprintf(2,'WARNING: Spectral integral frequency resolution estimation failed - defaulting to autocorrelation estimate');
			[fres,ierr] = var2fres(ARA,V);  % use autocorrelation-based estimate
		end
	end
	CAK = ARA;                              % CAK sequence for pre-optimisation
	H = var2trfun(ARA,fres);                % transfer function
	mdescript = sprintf('%d-variable VAR(%d)',n,r);
else
	[A0,C0,K0] = iss_rand(n,r,rho);         % random ISS model
	gc = ss_to_pwcgc(A0,C0,K0,V0);          % causal graph
	[A,C,K,V] = transform_ss(A0,C0,K0,V0);  % transform model to decorrelated-residuals form
	if isempty(fres)
		[fres,ierr] = ss2fres(A,C,K,V,[sitol,siminp2,simaxp2]);
		if isnan(fres) % failed!
			fprintf(2,'\nWARNING: Spectral integral frequency resolution estimation failed - defaulting to autocorrelation estimate');
			[fres,ierr] = ss2fres(A,C,K,V); % use autocorrelation-based estimate
		end
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
	gfile = fullfile(tempdir,'sim_model_pwcgc');
	wgraph2dot(n,eweight,gfile,[],gvprog,gvdisp);
	fprintf('\n');
end

% Save model

clearvars -except moddir fnamem mdescript V0 CAK H
wsfilem = fullfile(moddir,fnamem);
fprintf('*** saving model in ''%s''... ',wsfilem);
save(wsfilem);
fprintf('done\n\n');
