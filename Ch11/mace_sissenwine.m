function mace_sissenwine

% Re-create the SPR histogram from Mace & Sissenwine (1993)

SPR = 5:5:70;
Cum = [10 30 45 57 72 77 82 86 92 93 95 97 99 100];

figure(1)
clf
set(gcf,'units','cent','position',[10,10,9,8])

hold on
bar(SPR,Cum,1,'facecolor',repmat(0.8,[1,3]))
plot([0 100],[80 80],'k--')

set(gca,'tickdir','out','ticklength',[0.02 0.02],...
    'xlim',[0 75],'ylim',[0 100],'xtick',0:10:70,...
    'fontsize',12)
xlabel('Replacement SPR (%)','fontsize',14)
ylabel('Cumulative frequency (%)','fontsize',14)

