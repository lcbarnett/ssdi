function [CAK,H,fres,V0,mdescript] = sim_model(parms,mfile,gvparms,verb)

% Set up SS or VAR model, calculate transfer function, CAK sequence, etc.

if nargin < 2, || isempty(mfile), mfile = false; end
gvgc = nargin > 2 && ~isempty(gvparms);
if nargin < 4 || isempty(verb), verb = true; end

varmod = isfield(parms,'G'); % else SS model
if parms.varmod
	if isscalar(parms.G), parms.G = ones(parms.G); end;
	parms.n = size(G,1);
	assert(isfield(parms,'p'));
	assert(isfield(parms,'w'));
else
	assert(isfield(parms,'n'));
	assert(isfield(parms,'r'));
end
assert(isfield(parms,'fres'));
assert(isfield(parms,'rho' ));
assert(isfield(parms,'rmii'));

V0 = corr_rand(parms.n,parms.rmii); % residuals covariance matrix
if varmod
	A0 = var_rand(parms.G,parms.p,parms.rho,parms.w);
	if gvgc
		gc = var_to_pwcgc(A0,V0);
	end
	[A,V] = transform_var(A0,V0);
	if isempty(parms.fres)
		fres = var2fres(A);
	else
		fres = parms.fres;
	end
	CAK = A;
	H = var2trfun(A,fres);
	mdescript = sprintf('%d-variable VAR(%d)',parms.n,parms.p);
else
	[A0,C0,K0] = iss_rand(parms.n,parms.r,parms.rho);
	if gvgc
		gc = ss_to_pwcgc(A0,C0,K0,V0);
	end
	[A,C,K,V] = transform_ss(A0,C0,K0,V0);
	if isempty(parms.fres)
		fres = ss2fres(A,C,K);
	else
		fres = parms.fres;
	end
	CAK = iss2cak(A,C,K);
	H = ss2trfun(A,C,K,fres);
	mdescript = sprintf('%d-variable ISS(%d)',parms.n,parms.p);
end

if verb > 0
	fprintf('--------------------------------------------\n');
	fprintf('Model                : %s\n',mdescript);
	fprintf('--------------------------------------------\n');
	fprintf('Dimension            : %d\n',n);
	fprintf('Complexity (CAK)     : %d x %d x %d\n',size(CAK,1),size(CAK,2),size(CAK,3));
	fprintf('Frequency resolution : %d\n',fres);
	fprintf('--------------------------------------------\n');
end

if gvgc
	assert(isfield(gvparms,'file'));
	assert(isfield(gvparms,'prog'));
	assert(isfield(gvparms,'disp'));
	if isempty(gvparms.file), gvparms.file = fullfile(tempdir,'sim_model'); end
	eweight = gc/nanmax(gc(:));
	wgraph2dot(n,eweight,gvparms.file,[],gvparms.prog,gvparms.disp);
end

if ~isscalar(mfile) || mfile
	if isempty(mfile), mfile = fullfile(tempdir,'sim_model.mat'); end
	fprintf('*** saving model in ''%s''... ',mfile);
	save(mfile,'parms','CAK','H','fres','V0','mdescript');
	fprintf('done\n');
end
