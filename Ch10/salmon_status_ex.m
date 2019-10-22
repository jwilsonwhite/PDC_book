function salmon_status_ex

% Create example extinction risk figure, like in McClure et al. 2006

maxA = 8e3;
A = linspace(0,maxA,1e2); % abundance
P = linspace(0,12,1e2); % productivity

% Contours levels to plot
Con1 = 1*maxA./P;
Con2 = 3*maxA./P;
Con3 = 6*maxA./P;

% Create a bivariate function of P and A for stock status
X1 = (normpdf(log(P),log(4),0.3));
X2 = (normpdf(log(A),log(2.5e3),0.3));

XX = X1(:)*X2(:)';

% Plotting
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 8])

hold on
contourf(P,A,XX); % contour of stock status
plot(P,Con1,'k'); % contours of p(E)
plot(P,Con2,'k');
plot(P,Con3,'k');
colormap(1-gray)
set(gca,'ylim',[0 maxA])

set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',10)
xlabel('Productivity (recruits per spawner)','fontsize',12)
ylabel('Population abundance (no. spawners)','fontsize',12)

ch = colorbar;
set(ch,'ycolor','k','tickdir','out','ticklength',0.01,'fontsize',10)
ylabel(ch,'Probability density','fontsize',12)
