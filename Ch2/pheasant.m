function pheasant

% Plot growth of ring-neck pheasant population at Protection island. 
% Data from Einarsen (1945) Murrelet 26: 2-9

% The data:
Years = 1937.25:0.5:1942.75; % semiannual survey (March & November)
N = [8, 40, 30, 100, 81, 426, 282, 844, 705, 1540, 1325, 1898];

% Approximate model fit:
b = regress(log10(N(:)),[ones(length(N),1),Years(:)]);
n0 = b(1)+b(2)*1937.25;
r = 10^b(2);
Nm = (10^n0)*r.^(Years-Years(1));

% Setup figure
figure(1)
set(gcf,'units','cent','position',[40,40,8,8])
clf

% Plot on original axes
axes('position',[0.25 0.13 0.7 0.775])
hold on
plot(Years,Nm,'linestyle','-','color',[0.6 0.6 0.6])
plot(Years,N,'linestyle','-','color','k')
set(gca,'ylim',[0 2000],'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ytick',0:500:2000,'fontsize',10)
set(gca,'xlim',[1937,1943.5],'xtick',1937:1:1943,'xticklabel',{'1937','','1939','','1941','','1943'})
ylabel(gca,'Population size (N_t)','fontsize',12)

% Plot on semilog axes in inset
axes('position',[0.36 0.6 0.26 0.26])
hold on
plot(Years,log10(Nm),'-','color',[0.6 0.6 0.6])
plot(Years,log10(N),'linestyle','-','color','k')
set(gca,'ylim',[0.5 3.5],'tickdir','out','ticklength',[0.03 0.03])
set(gca,'ytick',0:1:4,'fontsize',8)
set(gca,'xlim',[1937,1943.5],'xtick',1937:1:1943,'xticklabel',{'1937','','','1940','','','1943'})
ylabel(gca,'log_1_0(N_t)','fontsize',8)




