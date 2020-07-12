% InfoVals
% Information values - plot to determine the Fisher information for various
% combinations of <V> and \hat{D}
% Assume constant error for these sets of plots
% Author: PKR, UNM, November 2015

close all;
clear;

% seed the random numbers so this code always produces the exact same
% results
rng(8675309);

SE = 1; % standard error of observations (uniformly variable distributions)
steps = 21;

T = (1:steps)';
medpath = SE*ones(steps,1);

Dsamp = logspace(log10(1e-4),log10(1e4),500)';
    
% get DST Information without motion blur
conFisher = recursiveFisher1D(Dsamp, T, medpath, 0);

% get DST Information with full frame motion blur
conFisherb = recursiveFisher1D(Dsamp, T, medpath, 1);

% plot the informations
figure;
plot(log(Dsamp),Dsamp.^2.*conFisher,'k','LineWidth',3);
hold on
plot(log(Dsamp),Dsamp.^2.*conFisherb,'k--','LineWidth',3);
plot([log(Dsamp(1)) log(Dsamp(end))],[(steps-1)/2 (steps-1)/2],'k--','LineWidth',1);
ax = gca;
fh = gcf;
set(ax,'FontSize',14,'FontWeight','bold');
xlabel('\boldmath$\hat{D}t/ \langle V \rangle$','FontSize',14,'Interpreter','latex')
ylabel('\boldmath$\tilde{\mathcal{I}}(\ln \hat{D}) $','FontSize',14,'Interpreter','latex')
legend({'\boldmath $t_{\epsilon} = 0$ ','\boldmath $t_{\epsilon} = \delta t$'},...
    'Interpreter','latex','Location','northwest')
ax.YTick = [0 (steps-1)/2];
ax.YTickLabel = {'0','(N-1)/2'};
ax.XTick = [-2 4];
ax.XTickLabel = {'\alpha','\beta'};
plot([-2 -2],[0 steps/2],'k--','LineWidth',1);
plot([4 4],[0 steps/2],'k--','LineWidth',1);
ax.YLim = [0 steps/2];
ax.XLim = [-9 9];
hold off

%% NOTE: I had to manually move the y-label because it auto sets outside of the tick label!

% scale for pdf
set(fh,'Units','Inches');
pos = get(fh,'Position');
set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh,'qualitate_info','fig');
print(fh,'qualitate_info','-dpdf','-r0')
print(fh,'qualitate_info','-dpng');


%% plot various trajectory lengths for constant SE
Len = [2 3 10 100];
cpath = cell(length(Len),1);
fishplot = cell(length(Len),1);
lfishplot = cell(length(Len),1);
for ii = 1:length(Len)
    tempT = (1:Len(ii))';
    cpath{ii} = SE*ones(Len(ii),1);
    fishplot{ii} = recursiveFisher1D(Dsamp, tempT, cpath{ii}, 1);
    lfishplot{ii} = fishplot{ii}.*Dsamp.^2;
end

% plot the various information curves for various track lengths
figure;
hold on;
dlsamp = log(Dsamp);
plot(dlsamp, lfishplot{1}/(Len(1)-1)*2,'k','LineWidth',2)
plot(dlsamp, lfishplot{2}/(Len(2)-1)*2,'k--','LineWidth',2)
plot(dlsamp, lfishplot{3}/(Len(3)-1)*2,'k:','LineWidth',2)
plot(dlsamp, lfishplot{4}/(Len(4)-1)*2,'k-.','LineWidth',2)
qq = strtrim(cellstr(num2str(Len'))');
legstring = strcat({'N = '},qq);
legend(legstring,'Location','northwest')

plot([dlsamp(1), dlsamp(end)],[1 1],'k--','LineWidth',1)
ax = gca;
fh = gcf;
set(ax,'FontSize',14,'FontWeight','bold');
xlabel('\boldmath$\hat{D}t/ \langle V \rangle$','FontSize',14,'Interpreter','latex')
ylabel('\boldmath$\tilde{\mathcal{I}}(\ln \hat{D}) $','FontSize',14,'Interpreter','latex')

ax.YTick = [0 1];
ax.YTickLabel = {'0','(N-1)/2'};
ax.YLim = [0 1.1];
ax.XLim = [-9 9];
ax.XTick = [-4.6052 0 4.6052];
ax.XTickLabel = {'10^{-2}', '10^0', '10^2'};
hold off;

%% NOTE: I had to manually move the y-label because it auto sets outside of the tick label!

% scale for pdf
set(fh,'Units','Inches');
pos = get(fh,'Position');
set(fh,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(fh,'quantitate_mean_info','fig');
print(fh,'quantitate_mean_info','-dpdf','-r0')
print(fh,'quantitate_mean_info','-dpng');