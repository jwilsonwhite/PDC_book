function SRcurves_plot

% Plot example of S-R relationships to illustrate location of cohort
% resonance & 2T cycle behavior

S = linspace(0,1,1e2);
a = 0.5;
b = 10;
aR = 0.5;
bR = 5;
Rick = aR.*S.*exp(-bR.*S);
Bev = a.*S./(1 + b.*S);

figure(1)
clf
set(gcf,'units','cent','position',[10 10 9,18])
sh(1) = subplot(2,1,1);
hold on
plot(S,Bev.*b/a,'k-')
plot([0 1],[0 1/0.05],'k--')
plot([0 1],[0 1/0.5],'k--')
plot([0 1],[0 1],'k--')
ylim([0 1])


sh(2) = subplot(2,1,2);
hold on
plot(S,Rick/(aR/bR*exp(-1)),'k-')
plot([0 1],[0 1/0.05],'k--')
plot([0 1],[0 1/0.5],'k--')
plot([0 1],[0 1],'k--')
ylim([0 1])

set(sh,'tickdir','out','ticklength',[0.015 0.015])
set(sh,'xtick',0:0.2:1,'ytick',0:0.2:1)
set(sh,'fontsize',10)
xlabel(sh(2),{'Egg production';'(proportion of unharvested maximum'},'fontsize',14)
ylabel(sh(1),{'Recruit abundance';'(proportion of maximum)'},'fontsize',14)

