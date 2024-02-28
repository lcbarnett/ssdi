%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defvar('n',        20   ); % observables dimension
defvar('r',        7    ); % state dimension (model order)
defvar('rho',      0.9  ); % VAR spectral radius
defvar('alpha',    0.05 ); % significance levels
defvar('mhtc',     true ); % apply multiple hypothesis test correction?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: need r < n for generic white-noise and causal cores to exist

[~,C,K] = iss_rand(n,r,rho);   % random ISS model (only need C and K parameters

[LC,LK] = iss_perfect_dd(C,K); % LC = white-noise core, LK =  causal cores

if mhtc, nhyp = 2*n; else, nhyp = 1; end % number of hypothesis: n x (left-tail and right-tail)
slev(1) = alpha/nhyp;   % Bonferroni correction for 2*n hypotheses (2 per channel) - left  tail
slev(2) = 1-alpha/nhyp; % Bonferroni correction for 2*n hypotheses (2 per channel) - right tail

thetaC = (pi/2)*gmetrics1(LC); % angles (in radians) of channel axes with white-noise core hyperplane
thetaK = (pi/2)*gmetrics1(LK); % angles (in radians) of channel axes with causal core hyperplane

cvals = get_haxa_cvals(n,[n-r,r],slev); % NOTE: dim(LC) = n-r, dim(LK) = r

% Compare with critical values - left tail

hC1 = thetaC < cvals(1,1); % significant channel contribution (i.e., reject H0) for white-noise core
hK1 = thetaK < cvals(2,1); % significant channel contribution (i.e., reject H0) for causal core

% Compare with critical values - right tail

hC2 = thetaC > cvals(1,2); % significant channel non-contribution (i.e., reject H0) for white-noise core
hK2 = thetaK > cvals(2,2); % significant channel non-contribution (i.e., reject H0) for causal core

fprintf('\nSignificantly small (left tail)\n   C   K\n-----------\n');
disp([hC1 hK1]);

fprintf('\nSignificantly large (right tail)\n   C   K\n-----------\n');
disp([hC2 hK2]);

% Display with confidence regions

aticks = [0 pi/4 pi/2];
aticklabels = {'0','\pi/4','\pi/2'};
xlims = [0.3,n+0.7];

figure(1); clf;

% White-noise core

subplot(2,1,1);
bar(thetaC);
patch([0 n+1 n+1 0], [cvals(1,1) cvals(1,1) cvals(1,2) cvals(1,2)],'k','FaceAlpha',0.1);
yline(cvals(1,1),'r');
yline(cvals(1,2),'g');
title('White-noise core axis angles');
xlim(xlims);
xlabel('channel');
xticks(1:n);
xticklabels(1:n);
ylim([0 pi/2]);
ylabel('radians');
yticks(aticks);
yticklabels(aticklabels);

% Causal core

subplot(2,1,2);
bar(thetaK);
patch([0 n+1 n+1 0], [cvals(2,1) cvals(2,1) cvals(2,2) cvals(2,2)],'k','FaceAlpha',0.1);
yline(cvals(2,1),'r');
yline(cvals(2,2),'g');
title('Causal core axis angles');
xlim(xlims);
xlabel('channel');
xticks(1:n);
xticklabels(1:n);
ylim([0 pi/2]);
ylabel('radians');
yticks(aticks);
yticklabels(aticklabels);
