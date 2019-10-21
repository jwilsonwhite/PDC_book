function Gaines_barnacles

% Replot data from Gaines&Roughgarden 1985


% Load data:
KLM_S = importdata('Gaines_data/KLM_S.csv');
KLM_F = importdata('Gaines_data/KLM_freespace.csv');
PR_S = importdata('Gaines_data/PR_S.csv');
PR_F = importdata('Gaines_data/PR_F.csv');

% Panel 1: KLM

% Settlement:
% Last two datapoints are the bottom of the x-axis
Date = KLM_S.data(19:end,6);
S = KLM_S.data(19:end,7);
S0 = S(end);
S = S0-S;
D0 = Date(end);
Dmax = Date(end-1);
Date = (D0-Date)/(D0-Dmax)*103; % convert to weeks
DateS = Date;

% Freespace:
Date = KLM_F.data(:,6);
F = KLM_F.data(:,7);
F0 = F(end);
F = F0-F;
F = F/100;
D0 = Date(end);
Dmax = Date(end-1);
Date = (D0-Date)/(D0-Dmax)*103; % convert to weeks
DateF = Date;

figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 18])
subplot(2,1,1)
hold on
plot(DateF(1:end-2),F(1:end-2),'k','linestyle','-')
set(gca,'xtick',0:20:100,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ytick',0:0.2:1)
set(gca,'xlim',[0 103],'ylim',[0 1],'fontsize',10)
Pos =  get(gca,'position');
axes('position',Pos)
plot(DateS(1:end-2),S(1:end-2),'k','linestyle','--')
set(gca,'color','none','yaxislocation','right','fontsize',10)
set(gca,'xtick',[],'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ytick',0:0.2:0.8)
set(gca,'xlim',[0 103],'ylim',[0 0.8])

% Panel 2: Pete's Rock

% Settlement:
% Last two datapoints are the bottom of the x-axis
Date = PR_S.data(1:end,6);
S = PR_S.data(1:end,7);
S0 = S(end);
S = S0-S;
S = max(S,0);
D0 = Date(1);
Dmax = Date(end-1);
Date = (D0-Date)/(D0-Dmax)*103; % convert to weeks
DateS = Date;

% Freespace:
Date = PR_F.data(:,6);
F = PR_F.data(:,7);
F0 = F(end);
F = F0-F;
D0 = Date(end);
Dmax = Date(end-1);
Date = (D0-Date)/(D0-Dmax)*103; % convert to weeks
DateF = Date;

subplot(2,1,2)
hold on
plot(DateF(1:end-2),F(1:end-2),'k','linestyle','-')
set(gca,'xtick',0:20:100,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ytick',0:0.2:1)
set(gca,'xlim',[0 103],'ylim',[0 1],'fontsize',10)
Pos =  get(gca,'position');
xlabel('Time (weeks)','fontsize',12)
%ylabel('Free space (proportion of total area','fontsize',12)
axes('position',Pos)
plot(DateS(1:end-2),S(1:end-2),'k','linestyle','--')
set(gca,'color','none','yaxislocation','right','fontsize',10)
set(gca,'xtick',[],'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ytick',0:3:15)
set(gca,'xlim',[0 103],'ylim',[0 15])
