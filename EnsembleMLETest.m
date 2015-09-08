% EnsembleMLETest
% Script to test the DEstimator ensemble methods

% Generate a cell array of trajectories of random sizes from 5 to 30
N = 30; % number of ensemble Trajectories
Nlen = round(5 + 25*rand(N,1)); % vector of trajectory lengths for each cell array

% Determine other simulation parameters
Dsim = 1;
dT = 1;
meanSE = 0.5; % use a constant standard error for this test
exposureT = dT; % full frame integration
alpha = exposureT/2/dT;
blurVar = 2*Dsim*exposureT*(1/3-alpha/2); % variance of frame averaged position due to full frame continuous exposure motion blur

Trajectories = cell(N,1);
% loop over each trajectory and generate it given the trajectory length
for ii = 1:N
   xpos = sqrt(2*Dsim*dT)*randn(Nlen(ii)+1,1);
   xpos = cumsum(xpos); % get true positions
   xmean = (xpos(1:end-1)+xpos(2:end))/2; % get expected averages
   % find variance due to motion blur and localization to get observations
   Obs = xmean + sqrt(meanSE^2+blurVar)*randn(Nlen(ii),1);
   T = dT*(0:Nlen(ii)-1)'; % get start frame times
   SE = meanSE*ones(Nlen(ii),1); % using constant errors for this sim
   Trajectories{ii} = [Obs T SE];
end

% Figure out the MLE positions
[D_mle, D_mle_llh] = DEstimator.computeEnsembleMLE(Trajectories, exposureT);

% create a sample vector of D values in log space
DVec = logspace(log10(0.1),log10(10),200);

% calculate the log likelihood of the sampled D values
llh_vec = DEstimator.computeEnsembleLLH(DVec,Trajectories,exposureT);

% convert the log likelihood into a re-scaled likelihood
lh_vec = exp(llh_vec - max(llh_vec));

%Plot the results.
plot(DVec,lh_vec,'r-') %plot likelihood curve
hold('on');
plot(D_mle, exp(D_mle_llh-max(llh_vec)),'bs','MarkerFaceColor','k') %plot D_mle value location
title('Diffusion constant likelihood');
xlabel('Diffusion constant (um^2/s)');
ylabel('Relative Likelihood');
set(gca,'XScale','log')
legend('Relative Likelihood of D-value',...
    sprintf('Maximum likelihood estimate: %f (um^2/s)',D_mle));
