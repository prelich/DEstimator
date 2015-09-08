% MLE_X_Test
% Script to test the validity of returning the maximum likelihood positions
% with the Laplace method.
% Peter K. Relich (physx.grad@gmail.com)
clear;
close all; % clean up the workspace

% First set up the simulation parameters
Dsim = 1; % Simulated diffusion parameter
dT = 1; % Simulated trajectory frame time
N = 100; % Numer of N observations per track
meanSE = 0.5; % The mean variance term, constant variance for now
exposureT = dT; % Make the exposure time the same as frame time for now

% Generate a 1D simulated trajectory of true positions
X = DEstimator.simulate1DDiffusion(N+1, Dsim, dT); %make N+1 X values
% Convert true positions into observations
var_gen_method = 'constant'; % lets keep the localization variance constant for now

% I need to alter a few lines inplicit in the class function
% simulate1DObservationVar for testing purposes, hence the script is
% re-stated explicitly here rather invoked from the class
% Mainly the observation is going to be perturbed by random variables
% twice, which is the same as perturbing with random variables once with a
% larger standard deviation
if isscalar(meanSE)
    %If a scalar V is given
    SEvec = DEstimator.generateObsStandardError(meanSE,N,var_gen_method);
else
    SEvec = SEvec(1:N);
    SEvec = SEvec(:);
end
Tvec = (0:N-1)'.*dT;
alpha = exposureT/(2*dT);
Xbar_mean = (1-alpha)*X(1:end-1) + alpha*X(2:end);
Xbar_var = 2*Dsim*exposureT*(1/3-alpha/2);
Xbar_realize = Xbar_mean+randn(N,1).*sqrt(Xbar_var); % we want to know the realizations
ObSim = Xbar_realize + randn(N,1).*SEvec;

% [ObSim, Tvec, Vvec] = DEstimator.simulate1DObservationVar(X, meanVar, Dsim, dT, exposureT,var_gen_method);

% Create an instance of the DEstimator class
est = DEstimator(ObSim, Tvec, SEvec, dT);

% Perform the MLE estimation of the true positions
[MLE_X, MLE_SE, MLE_LLH_X] = est.MLEPositions(Dsim);

% Perform the MLE estimation of true positions with the direct method
[~, posMLE]= DEstimator.LLH_laplaceDirect1D(Dsim, ObSim, Tvec, SEvec, exposureT);

% Perform MLE estimation of D
[MLE_D, MLE_LLH_D] = est.MLE;

% Perform the MLE estimation of the true positions with the mle D
[MLE_XP, MLE_LLH_XP] = est.MLEPositions(MLE_D);

% Calculate the squared error lost on the MLE positions vs. actual
% positions
sqLoss_trueD = (MLE_X - X(1:end-1)).^2;
sqLoss_mleD = (MLE_XP - X(1:end-1)).^2;

sqLoss_trueD2 = (MLE_X - X(2:end)).^2;
sqLoss_mleD2 = (MLE_XP - X(2:end)).^2;

sqLoss_trueD3 = (MLE_X - Xbar_mean).^2;
sqLoss_mleD3 = (MLE_XP - Xbar_mean).^2;

sqLoss_trueD4 = (MLE_X - Xbar_realize).^2;
sqLoss_mleD4 = (MLE_XP - Xbar_realize).^2;

meanLossTrueD = mean(sqLoss_trueD)
meanLossmleD = mean(sqLoss_mleD)
meanLossTrueD2 = mean(sqLoss_trueD2)
meanLossmleD2 = mean(sqLoss_mleD2)
meanLossTrueD3 = mean(sqLoss_trueD3)
meanLossmleD3 = mean(sqLoss_mleD3)
meanLossTrueD4 = mean(sqLoss_trueD4)
meanLossmleD4 = mean(sqLoss_mleD4)

% Plot the various squared loss methods
figure;
plot(1:length(MLE_X),sqLoss_trueD,'LineWidth',2,'DisplayName',...
    'Loss from frame start positions');
hold on
% plot(1:length(MLE_X),sqLoss_mleD,'--','LineWidth',2);
plot(1:length(MLE_X),sqLoss_trueD2,'LineWidth',2,'DisplayName',...
    'Loss from frame end positions');
% plot(1:length(MLE_X),sqLoss_mleD2,'--','LineWidth',2);
plot(1:length(MLE_X),sqLoss_trueD3,'LineWidth',2,'DisplayName',...
    'Loss from expected frame average positions');
% plot(1:length(MLE_X),sqLoss_mleD3,'--','LineWidth',2);
plot(1:length(MLE_X),sqLoss_trueD4,'LineWidth',2,'DisplayName',...
    'Loss from realized frame average positions');
hold off
lh = legend('show','location', 'northwest');
title('Squared Error Loss for MLE of various particle positions with Laplace','FontSize',14);
xlabel('Associated Trajectory Observation','FontSize',14);
ylabel('Squared Error Loss','FontSize',14);

%% Perform this same test with the laplace direct method!
% Calculate the squared error lost on the MLE positions vs. actual
% positions
DsqLoss = (posMLE - X).^2;

meanLossmleDF = mean(DsqLoss)

% Plot the various squared loss methods
figure;
plot(1:length(posMLE),DsqLoss,'LineWidth',2,'DisplayName',...
    'Loss from frame start positions');
hold on

lh = legend('show','location', 'northwest');
title('Squared Error Loss for MLE of various particle positions with Direct Laplace','FontSize',14);
xlabel('Associated Trajectory Observation','FontSize',14);
ylabel('Squared Error Loss','FontSize',14);