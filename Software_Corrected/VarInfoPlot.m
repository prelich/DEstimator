% VarInfoPlot
% script to plot the fisher information for variable localization errors
% Author: PKR, UNM, November 2015

close all;
clear;

% seed the random numbers so this code always produces the exact same
% results
rng(8675309);

SE = 1; % standard error of observations (uniform variable distributions)
samples = 1e3; % number of trajectory samples to produce
steps = 11; % number of observations in a trajectory

T = (1:steps)';
medpath = SE*ones(steps,1);
Dsamp = logspace(log10(1e-4),log10(1e4),500)';

% get constant error Information with full frame motion blur
conFisher = Dsamp.^2.*recursiveFisher1D(Dsamp, T, medpath, 1);

% Initialize vector on
rFisher = zeros(length(Dsamp),samples);

bound = 0.9; % normal distribution goes +- 80% of the mean

for ii = 1:samples
    Serr = (1 + bound*(2*rand(steps,1)-1));
    % make sure the mean of Serr^2 is SE^2
    Serr = SE*sqrt(Serr.^2/mean(Serr.^2));
    
    % Get recursive Information
    rFisher(:,ii) = Dsamp.^2.*recursiveFisher1D(Dsamp, T, Serr, 1);
end

figure;
plot(log(Dsamp),max(rFisher,[],2),'k--','LineWidth',3);
hold on
plot(log(Dsamp),conFisher,'k','LineWidth',2);
plot(log(Dsamp),min(rFisher,[],2),'k--','LineWidth',3);
legend({'Information Bounds','\langleV\rangle Information'},'Location','northwest');
ax = gca;
fh = gcf;
set(ax,'FontSize',14,'FontWeight','bold');
xlabel('\boldmath$\hat{D}t/ \langle V \rangle$','FontSize',14,'Interpreter','latex')
ylabel('\boldmath$\tilde{\mathcal{I}}(\ln \hat{D}) $','FontSize',14,'Interpreter','latex')
ax.YTick = [0 (steps-1)/2];
ax.YTickLabel = {'0','(N-1)/2'};
ax.XTick = 0;
ax.XTickLabel = '1';
ax.YLim = [0 steps/2];
ax.XLim = [-9 9];

%% NOTE: I had to manually move the y-label because it auto sets outside of the tick label!

% scale for pdf
set(fh,'Units','Inches');
pos = get(fh,'Position');
set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh,'extrema_info','fig');
print(fh,'extrema_info','-dpdf','-r0')
print(fh,'extrema_info','-dpng');
