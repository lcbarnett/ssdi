%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('n',        20    ); % observables dimension
defvar('r',        7     ); % state dimension (model order)
defvar('rho',      0.9   ); % VAR spectral radius
defvar('alpha',    0.05  ); % significance levels
defvar('mhtc',     true  ); % use (Bonferroni) multiple hypothesis test correction?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: need r < n for generic white-noise and causal cores to exist

[~,C,K] = iss_rand(n,r,rho);       % random ISS model (only need C and K parameters

[L_wn,L_cc] = iss_perfect_dd(C,K); %  L_cc =  causal core, L_wn = white-noise core

% Calculate Beta statistics

beta_cc = habeta(L_cc); % hyperplane/channel Beta statistics for causal core
beta_wn = habeta(L_wn); % hyperplane/channel Beta statistics for white-noise core

% Statistical inference on Beta statistics for all channels
%
% NOTE: dim(L_wn) = n-r, dim(L_cc) = r

[cval_cc,pval_cc,sig_cc] = habeta_statinf(beta_cc,n,r,  alpha,'both',mhtc);
[cval_wn,pval_wn,sig_wn] = habeta_statinf(beta_wn,n,n-r,alpha,'both',mhtc);

fprintf('\nSignificant non-participation (left tail):\n\n     CC    WN\n     --------\n');
disp(0+[sig_cc(:,1) sig_wn(:,1)]);

fprintf('Significant participation (right tail):\n\n     CC    WN\n     --------\n');
disp(0+[sig_cc(:,2) sig_wn(:,2)]);

% Display with confidence regions

figure(1); clf;
sgtitle(sprintf('$\\beta$ statistics for causal and white-noise cores: n = %d, r = %d\n',r,n),'Interpreter','latex','FontSize',16);
xlims = [0.3,n+0.7];

% Causal core

subplot(2,1,1);
bar(beta_cc);
patch([0 n+1 n+1 0], [cval_cc(1) cval_cc(1) cval_cc(2) cval_cc(2)],'k','FaceAlpha',0.1);
yline(cval_cc(1),'r','LineWidth',1.2);
yline(cval_cc(2),'r','LineWidth',1.2);
title('Causal core');
xlim(xlims);
xlabel('channel');
xticks(1:n);
xticklabels(1:n);
ylim([0 1]);
yticks(0:0.2:1);
ylabel('\beta statistic');

% White-noise core

subplot(2,1,2);
bar(beta_wn);
patch([0 n+1 n+1 0], [cval_wn(1) cval_wn(1) cval_wn(2) cval_wn(2)],'k','FaceAlpha',0.1);
yline(cval_wn(1),'r','LineWidth',1.2);
yline(cval_wn(2),'r','LineWidth',1.2);
title('White-noise core');
xlim(xlims);
xlabel('channel');
xticks(1:n);
xticklabels(1:n);
ylim([0 1]);
yticks(0:0.2:1);
ylabel('\beta statistic');
