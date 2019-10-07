function chinook

% Plot population trajectory of Sacramento R. Chinook salmon.
% Data digitized using ImageJ from Botsford & Brittnacher (1998) Con Bio

% The data:
D = importdata('Botsford_Brittnacher_chinook_data.csv');
Y = D.data(:,1);
N = D.data(:,2);


% Approximate model fit by regression on log(N)
b = regress(log10(N(:)),[ones(length(N),1), Y(:)]);
n0 = b(1)+b(2)*1968;
r = 10^(b(2));
Nm = (10^n0)*r.^(Y-1968);
%keyboard

% Setup figure
figure(2)
set(gcf,'units','cent','position',[40,40,8,8])
clf

% Plot on original axes
axes('position',[0.25 0.13 0.7 0.775])
hold on
plot(Y,Nm,'linestyle','-','color',[0.6 0.6 0.6])
plot(Y,N,'linestyle','-','color','k')
set(gca,'ylim',[0 1.2e5],'tickdir','out','ticklength',[0.02 0.02])
set(gca,'ytick',0:2e4:1.2e5,'fontsize',10)
set(gca,'xlim',[1965,1995],'xtick',1965:5:1995)
ylabel(gca,'Population size (N_t)','fontsize',12)

% Plot on semilog axes in inset
axes('position',[0.6 0.6 0.26 0.26])
hold on
plot(Y,log10(Nm),'-','color',[0.6 0.6 0.6])
plot(Y,log10(N),'linestyle','-','color','k')
set(gca,'ylim',[2 6],'tickdir','out','ticklength',[0.03 0.03])
set(gca,'ytick',0:1:6,'fontsize',8)
set(gca,'xlim',[1965,1995],'xtick',1965:5:1995)
ylabel(gca,'log_1_0(N_t)','fontsize',8)




