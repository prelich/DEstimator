% Trajectory_Sim
% Script to simulation trajectories and run diffusion estimates!
% WARNING: SCRIPT Takes ~8 hours to run on ~2012 Desktop PC architecture
% Author: PKR, UNM, Dec 2015

close all;
clear;
% make this script deterministic
rng(8675309);


%% Create an instance of the Simulation Class and Generate Trajectories
A = DSim; % Call the constructor

% Determine number of samples per batch  and the number of batches.  (DSim
% will run out of memory otherwise!)
samples = 1e2;
batches = 1e2;

% Determine simulation parameters
A.FrameSize = 100; % maximum frame length for a trajectory
A.MeanPhot = 200; % mean emission rate of photons per frame 
A.PsfSigma = 1; % point spread function sigma in pixels
A.Texp = 1; % exposure time
A.StateParams = [0.2 0.05]; % Kon and Koff states for a `blinking' probe

% Background parameters
A.BGsigma = 5;
A.BGTotalI = 2000;
A.Periodic = 21;

% Number of Trajectories to Sample
A.N = samples;

% figure out the conversion parameters from simulation to experiment
pixelSize = 0.1; % 0.1 um/pixel
frameTime = 0.01; % 100 frame/sec  
% Loop over orders of magnitude changes in D
Dsim = logspace(log10(1e-5),log10(1e-1),20);
% initialize analysis variables
pvals = [0.32 0.05]; % going to perform various confidence intervals
% for analysis that recognizes variable localization errors
MLEvar = zeros(A.N*batches,length(Dsim));
MLEvar_LLH = zeros(A.N*batches,length(Dsim));
obsInfvar = zeros(A.N*batches,length(Dsim));
logSigvar = zeros(A.N*batches,length(Dsim));
% for analysis that only recognizes a reduced localization error
MLEred = zeros(A.N*batches,length(Dsim));
MLEred_LLH = zeros(A.N*batches,length(Dsim));
obsInfred = zeros(A.N*batches,length(Dsim));
logSigred = zeros(A.N*batches,length(Dsim));
% for analysis that recognizes a reduced localization error after filtering
% out really poor fits
MLEfred = zeros(A.N*batches,length(Dsim));
MLEfred_LLH = zeros(A.N*batches,length(Dsim));
obsFInfred = zeros(A.N*batches,length(Dsim));
logSigfred = zeros(A.N*batches,length(Dsim));

Successvar = zeros(length(pvals),length(Dsim));
Successred = zeros(length(pvals),length(Dsim));
Successfred = zeros(length(pvals),length(Dsim));


for ll = 1:length(Dsim)
    A.D = Dsim(ll)*frameTime/pixelSize^2; % units of px^2/s
    % log truth
    lD = log(A.D);
    
    %% Perform D Estimate analysis    
    % loop over each trajectory
    for kk = 1:batches
        % Run the simulation and Return the trajectories for each batch
        A.runSim;
        [filterObs, filterVar] = A.filterCRLB;
        for ii = 1:A.N
            Obs = A.O{ii}(:,1:2);
            T = A.O{ii}(:,5);
            if isempty(T)
                continue;
            end
            SE = sqrt(A.CRLB{ii}(:,1:2));
            fObs = filterObs{ii}(:,1:2);
            fT = filterObs{ii}(:,5);
            % it has to be the mean of the variances to be more accurate for this type of analysis
            reducedSE = repmat(sqrt(mean(A.CRLB{ii}(:,1:2))),[length(T),1]);
            fSE = repmat(sqrt(mean(filterVar{ii}(:,1:2))),[length(fT),1]);
            % get maximum likelihoods
            index = (kk-1)*samples;
            zz = ii+index; % determine array value for each trajectory
            [MLEvar(zz,ll),MLEvar_LLH(zz,ll)]=computeMLE(Obs, T, SE, A.Texp);
            [MLEred(zz,ll),MLEred_LLH(zz,ll)]=computeMLE(Obs, T, reducedSE, A.Texp);
            [MLEfred(zz,ll),MLEfred_LLH(zz,ll)]=computeMLE(fObs, fT, fSE, A.Texp);
            % get observed err, information is additive, so add dimensions
            [~, varinftemp1]=LLH_ObsI_recursive1D(MLEvar(zz,ll), Obs(:,1), T, SE(:,1), A.Texp);
            [~, varinftemp2]=LLH_ObsI_recursive1D(MLEvar(zz,ll), Obs(:,2), T, SE(:,2), A.Texp);
            obsInfvar(zz,ll) = varinftemp1+varinftemp2;
            % additive information for the constant error case
            [~, coninftemp1]=LLH_ObsI_recursive1D(MLEred(zz,ll), Obs(:,1), T, reducedSE(:,1), A.Texp);
            [~, coninftemp2]=LLH_ObsI_recursive1D(MLEred(zz,ll), Obs(:,2), T, reducedSE(:,2), A.Texp);
            obsInfred(zz,ll) = coninftemp1 + coninftemp2;
            % additive information for the constant error, filtered case
            [~, coninfftemp1]=LLH_ObsI_recursive1D(MLEfred(zz,ll), fObs(:,1), fT, fSE(:,1), A.Texp);
            [~, coninfftemp2]=LLH_ObsI_recursive1D(MLEfred(zz,ll), fObs(:,2), fT, fSE(:,2), A.Texp);
            obsFInfred(zz,ll) = coninfftemp1 + coninfftemp2;
            
            for jj = 1:length(pvals)
                logSigvar(zz,ll)=waldInterval(MLEvar(zz,ll)^2*obsInfvar(zz,ll),pvals(jj));
                logSigred(zz,ll)=waldInterval(MLEred(zz,ll)^2*obsInfred(zz,ll),pvals(jj));
                logSigfred(zz,ll)=waldInterval(MLEfred(zz,ll)^2*obsFInfred(zz,ll),pvals(jj));
                % success rates for variable localization errors
                if isreal(logSigvar(zz,ll)) && lD < log(MLEvar(zz,ll)) ...
                        + logSigvar(zz,ll) && lD > log(MLEvar(zz,ll)) - logSigvar(zz,ll)
                    Successvar(jj,ll) = Successvar(jj,ll) + 1;
                end
                % success rates for constant localization errors
                if isreal(logSigred(zz,ll)) &&  lD < log(MLEred(zz,ll)) ...
                        + logSigred(zz,ll) && lD > log(MLEred(zz,ll)) - logSigred(zz,ll)
                    Successred(jj,ll) = Successred(jj,ll) + 1;
                end
                % success rates for constant localization errors with filtering
                if isreal(logSigfred(zz,ll)) &&  lD < log(MLEfred(zz,ll)) ...
                        + logSigfred(zz,ll) && lD > log(MLEfred(zz,ll)) - logSigfred(zz,ll)
                    Successfred(jj,ll) = Successfred(jj,ll) + 1;
                end
            end
        end
    end

end

% Get an estimate on the average CRLB value of all localizations from a
CRLBsample = cell2mat(A.CRLB);
CRLBsample = CRLBsample(:,1:2); % only want localiztions
CRLBsample = CRLBsample(:); % verticle vector
CRLBsample = CRLBsample(CRLBsample < 0.5); % remove huge errors
meanCRLB = mean(CRLBsample);

critfailvar = 0*Dsim;
critfailfred = 0*Dsim;
critfailred = 0*Dsim;

% find critical fail rates
for ii = 1:length(Dsim)
   critfailvar(ii) = nnz(imag(logSigvar(:,ii)));
   critfailfred(ii) = nnz(imag(logSigfred(:,ii)));
   critfailred(ii) = nnz(imag(logSigred(:,ii)));
end

norm = samples*batches; % normalize all success rates
xaxis = Dsim/meanCRLB*frameTime/pixelSize^2;
%% Perform plotting for the paper!

% plot success rates of confidence intervals as a function of Dt/CRLB;
figure;
hold on;
plot(log10(xaxis), Successvar(1,:)./(norm-critfailvar),'k-','LineWidth',2);
plot(log10(xaxis), Successfred(1,:)./(norm-critfailfred),'k--','LineWidth',2);
plot(log10(xaxis), Successred(1,:)./(norm-critfailred),'k-.','LineWidth',2);

legend({'VEA Interval','FMEA Interval',...
    'MEA Interval'},'Location','southeast');

plot(log10(xaxis), Successvar(2,:)./(norm-critfailvar),'k-','LineWidth',2);
plot(log10(xaxis), Successfred(2,:)./(norm-critfailfred),'k--','LineWidth',2);
plot(log10(xaxis), Successred(2,:)./(norm-critfailred),'k-.','LineWidth',2);

plot(log10(xaxis), 0*Dsim+0.68, 'k--');
plot(log10(xaxis), 0*Dsim+0.95, 'k--');
xlabel('D\tau/\langleV\rangle')
ylabel('Interval Success Rate');
ax = gca;
fh = gcf;
set(ax,'FontSize',14,'FontWeight','bold');
ax.YTick = [0.68 0.95];
ax.YTickLabel = {'68%','95%'};
ax.XTick = [-3 -2 -1 0];
ax.XTickLabel = {'10^{-3}','10^{-2}','10^{-1}','10^{0}'};
ax.YLim = [0.4 1];
ax.XLim = [-3 0.8];
hold off;

% scale for pdf
set(fh,'Units','Inches');
pos = get(fh,'Position');
set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh,'confidence_success','fig');
print(fh,'confidence_success','-dpdf','-r0')
print(fh,'confidence_success','-dpng');


% plot the interval widths for non-imaginary intervals
% remove all the failed/complex intervals
lSigvar = cell(length(Dsim),1);
lSigred = cell(length(Dsim),1);
lSigfred = cell(length(Dsim),1);
for ii = 1:length(Dsim)
    tempvar = logSigvar(:,ii);
    lSigvar{ii} = tempvar(~imag(tempvar));
    tempvar = logSigred(:,ii);
    lSigred{ii} = tempvar(~imag(tempvar));
    tempvar = logSigfred(:,ii);
    lSigfred{ii} = tempvar(~imag(tempvar));
end
meanlSigvar = 0*Dsim;
meanlSigred = 0*Dsim;
meanlSigfred = 0*Dsim;
% get the mean log sim values
for ii = 1:length(Dsim)
   meanlSigvar(ii) = mean(lSigvar{ii});
   meanlSigred(ii) = mean(lSigred{ii});
   meanlSigfred(ii) = mean(lSigfred{ii});
end

figure;
hold on
plot(log10(xaxis),meanlSigvar,'k-','LineWidth',2);
plot(log10(xaxis),meanlSigfred,'k--','LineWidth',2);
plot(log10(xaxis(4:end)),meanlSigred(4:end),'k-.','LineWidth',2);
plot(log10(xaxis),Dsim*0+0.5,'k--');
legend({'VEA Interval','FMEA Interval',...
    'MEA Interval','Fiducial'},'Location','northeast');

xlabel('D\tau/\langleV\rangle')
ylabel('\bf{Mean Interval Range} \boldmath{$|$}\bf{ln}\boldmath{$(\hat{D}/D)|$}','interpreter','latex');
ax2 = gca;
fh2 = gcf;
set(ax2,'FontSize',14,'FontWeight','bold');
ax2.XTick = [-3 -2 -1 0];
ax2.XTickLabel = {'10^{-3}','10^{-2}','10^{-1}','10^{0}'};
ax2.XLim = [-3 0.8];
hold off;

% scale for pdf
set(fh2,'Units','Inches');
pos = get(fh2,'Position');
set(fh2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh2,'confidence_widths','fig');
print(fh2,'confidence_widths','-dpdf','-r0')
print(fh2,'confidence_widths','-dpng');

% plot the critical failures next!
figure;
hold on
plot(log10(xaxis),critfailvar/norm,'k-','LineWidth',2)
plot(log10(xaxis),critfailfred/norm,'k--','LineWidth',2)
plot(log10(xaxis),critfailred/norm,'k-.','LineWidth',2)

legend({'VEA','FMEA',...
    'MEA'},'Location','northeast');

xlabel('D\tau/\langleV\rangle')
ylabel('Critical Failure Rate');
ax3 = gca;
fh3 = gcf;
set(ax3,'FontSize',14,'FontWeight','bold');
ax3.XTick = [-3 -2 -1 0];
ax3.XTickLabel = {'10^{-3}','10^{-2}','10^{-1}','10^{0}'};
ax3.XLim = [-3 0.8];
hold off

% scale for pdf
set(fh3,'Units','Inches');
pos = get(fh3,'Position');
set(fh3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh3,'critical_fail','fig');
print(fh3,'critical_fail','-dpdf','-r0')
print(fh3,'critical_fail','-dpng');


% plot the log variances of the MLE values
% remove all the failed MLE values
cMLEvar = cell(length(Dsim),1);
cMLEred = cell(length(Dsim),1);
cMLEfred = cell(length(Dsim),1);
for ii = 1:length(Dsim)
    tempvar = logSigvar(:,ii);
    tempMLE = MLEvar(:,ii);
    cMLEvar{ii} = tempMLE(~imag(tempvar));
    tempvar = logSigred(:,ii);
    tempMLE = MLEred(:,ii);
    cMLEred{ii} = tempMLE(~imag(tempvar));
    tempvar = logSigfred(:,ii);
    tempMLE = MLEfred(:,ii);
    cMLEfred{ii} = tempMLE(~imag(tempvar));
end
meanMLEvar = 0*Dsim;
meanMLEred = 0*Dsim;
meanMLEfred = 0*Dsim;
varup = 0*Dsim;
vardown = 0*Dsim;
redup = 0*Dsim;
reddown = 0*Dsim;
fredup = 0*Dsim;
freddown = 0*Dsim;

% get means of various MLE values
% also get 16% and 84% quantiles for 1 sigma deviation
for ii = 1:length(Dsim)
   meanMLEvar(ii) = mean(cMLEvar{ii});
   meanMLEred(ii) = mean(cMLEred{ii});
   meanMLEfred(ii) = mean(cMLEfred{ii});
   
   varup(ii) = quantile(cMLEvar{ii},0.975);
   vardown(ii) = quantile(cMLEvar{ii},0.025);
      
   redup(ii) = quantile(cMLEred{ii},0.975);
   reddown(ii) = quantile(cMLEred{ii},0.025);
   
   fredup(ii) = quantile(cMLEfred{ii},0.975);
   freddown(ii) = quantile(cMLEfred{ii},0.025);
end

% calculate the CRLB for constant error for various Dsim values, assume no
% intermittency, best information with constant error!
conFisher = recursiveFisher1D(Dsim, (1:A.FrameSize)', sqrt(meanCRLB)*ones(A.FrameSize,1), 1);
fidCRLB = 1./conFisher';

% plot the mean MLE's, their STD's and the CRLB
figure;
loDsim = log(Dsim);
hold on;
plot(log10(xaxis),(meanMLEvar-Dsim)./Dsim,'k-','LineWidth',2);
plot(log10(xaxis),(meanMLEfred-Dsim)./Dsim,'k--','LineWidth',2);

% fill spacing
spacingx = [log10(xaxis), log10(xaxis(end:-1:1))];

% variable error fill region!
varfill = [(varup-Dsim)./Dsim, (vardown(end:-1:1)-Dsim(end:-1:1))./Dsim(end:-1:1)];
q3=patch(spacingx,varfill, 'r');
% Reduced error fill region!
fredfill = [(fredup-Dsim)./Dsim, (freddown(end:-1:1)-Dsim(end:-1:1))./Dsim(end:-1:1)];
q1=patch(spacingx,fredfill, 'g');
% CRLB fill region!
CRLBup = (1.96*sqrt(fidCRLB))./Dsim;
CRLBdown = (-1.96*sqrt(fidCRLB))./Dsim;
CRLBfill = [CRLBup, CRLBdown(end:-1:1)];
q2=patch(spacingx, CRLBfill, 'b');

set(q1,'facealpha',.1,'edgealpha',0);
set(q2,'facealpha',.1,'edgealpha',0);
set(q3,'facealpha',.1,'edgealpha',0);

legend({'\bf{VEA} \boldmath{$\langle\hat{D}\rangle$}',...
    '\bf{FMEA} \boldmath{$\langle\hat{D}\rangle$}',...
    '\bf{VEA} \boldmath{$95\%$} \bf{Empirical Interval}',...
    '\bf{FMEA} \boldmath{$95\%$} \bf{Empirical Interval}',...
    '\bf{CRLB} \boldmath{$95\%$} \bf{Analytical Interval}'},'Interpreter','latex');

plot(log10(xaxis),Dsim*0,'k--');

xlabel('D\tau/\langleV\rangle')
ylabel('\boldmath$(\hat{D}-D)/D$','interpreter','latex');

ax4 = gca;
fh4 = gcf;
ax4.YLim = [-4 10];
ax4.XLim = [-3 0.8];
set(ax4,'FontSize',14,'FontWeight','bold');
ax4.XTick = [-3 -2 -1 0];
ax4.XTickLabel = {'10^{-3}','10^{-2}','10^{-1}','10^{0}'};

hold off

% scale for pdf
set(fh4,'Units','Inches');
pos = get(fh4,'Position');
set(fh4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh4,'MLE_var','fig');
print(fh4,'MLE_var','-dpdf','-r0')
print(fh4,'MLE_var','-dpng');

