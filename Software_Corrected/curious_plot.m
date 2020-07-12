% curious_plot
% Plot to see if Stein Loss works for arbitrary distributions

% rng(8675309);

SE = 1; % standard error of observations (uniform variable distributions)
steps = 21; % number of observations in a trajectory

T = (1:steps)';
medpath = SE*ones(steps,1);
Dsamp = logspace(log10(1e-2),log10(1e2),500)';

Dtrue = 1;

% get constant error Information with full frame motion blur
conFisher = Dsamp.^2.*recursiveFisher1D(Dsamp, T, medpath, 0);

bound = 0.9; % normal distribution goes +- 80% of the mean

% only one sample
Serr = (1 + bound*(2*rand(steps,1)-1));
% make sure the mean of Serr^2 is SE^2
Serr = SE*sqrt(Serr.^2/mean(Serr.^2));

% Get recursive Information
rFisher = Dsamp.^2.*recursiveFisher1D(Dsamp, T, Serr, 0);

% generate observations
Obs = cumsum(sqrt(2*Dtrue)*randn(steps,1))+Serr.*randn(steps,1);

% get likelihood distribution
LLH = LLH_recursive1D(Dsamp, Obs, 1:length(Obs), Serr, 0);
% compute the MLE
[MLE_val, LLH_MLE] = computeMLE(Obs, 1:length(Obs), Serr, 0);
% compute the observed information
[~, observed_info] = LLH_ObsI_recursive1D(MLE_val, Obs, 1:length(Obs), Serr, 0);
observed_info = observed_info*MLE_val^2;

% Get recursive Information at point
MLErFisher = MLE_val^2.*recursiveFisher1D(MLE_val, T, Serr, 0);

% define the stein loss and quadratic approximation
logDrat = log(MLE_val./Dsamp);
quadratic = observed_info*logDrat.^2/2;
stein = observed_info*(exp(logDrat) - logDrat - 1);

% plot all functions
figure;
hold on;
plot(logDrat, LLH_MLE-LLH,'Linewidth',2);
plot(logDrat, quadratic, '--','Linewidth',2);
plot(logDrat, stein, '-.','LineWidth',2);
xlabel('Sampled D value','FontSize',16);
ylabel('LLR','Fontsize',16);
legend({'Actual LLH','Quadratic LLH','Stein LLH'});
hold off
