function urchin_growth_mortality_plot

% Plot the hypothetical growth & mortality functions used in urchin model

L = 0:100; % length
D = 0.60 - 0.59/50^2*(L-50).^2; % min at 0.01, peak of 0.6 at 50 mm 
dg = ones(length(L),1)*0.3; % absolute value of length


figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 6])
hold on

plot(L,D,'k-','linewidth',1)
plot(L,dg,'k--','linewidth',1)

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10)
ylim([0 0.7])
set(gca,'xtick',0:20:100,'xlim',[0 100])
ylabel('|dg/dt| and |D|','fontsize',12)
xlabel('Test diameter (mm)','fontsize',12)


