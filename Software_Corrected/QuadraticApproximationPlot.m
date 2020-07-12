% QuaddraticApproximationPlot
% Quadratic Approximation Limit for Diffusion estimation on No Error trajectories 

logDrat = -1.5:0.01:1.5;
response = exp(logDrat) - logDrat - logDrat.^2/2 - 1;

plot(logDrat,response,'k','LineWidth',2);
hold on;
plot(logDrat,0*logDrat,'--k','LineWidth',2);
f = gcf;
ha = gca;
set(ha,'FontSize',14,'FontWeight','bold');
xlabel('\bf{ln} \boldmath$\hat{D}/D$','Interpreter','latex');
ylabel('LLR / Information');
legend({'Deviation between Approximate and Exact LLRs','Fiducial for Perfect Agreement'},...
    'Location','northwest')
hold off;

% plot the quadratic super-imposed over the Stein loss
quadratic = logDrat.^2;
stein = 2*(exp(logDrat) - logDrat - 1);

Drat = exp(logDrat);
nappx = (Drat-1).^2;

figure;
plot(logDrat,stein,'k','LineWidth',2);
hold on;
plot(logDrat,quadratic,'k--','LineWidth',2);
plot(logDrat,nappx,'k-.','LineWidth',2);
f2 = gcf;
ha2 = gca;
set(ha2,'FontSize',14,'FontWeight','bold');
xlabel('\bf{ln} \boldmath$\left( \hat{D}/D \right)$','Interpreter','latex');
ylabel('LLR / Information');
legend({'Exact LLR Distribution','Log-normal Approximation','Normal Approximation'},'Location','northwest');
plot([-.5 -.5],[4 0],'k--');
plot([.5 .5],[4 0],'k--');
ha2.YLim = [0 4];
hold off;

% scale for pdf
set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

% save the plot
saveas(f2,'quad_approx_noerr','fig');
print(f2,'quad_approx_noerr','-dpdf','-r0')
print(f2,'quad_approx_noerr','-dpng');