function ricker_bevholt

% Plot Ricker & Beverton-Holt curves

a = [1.1 2 4]; % values of a parameter
b = 0.1; % b for Ricker
b2 = 10; % b for BH
N = 0:0.01:100; % initial N values
LS = {'-','--',':'}; % line styles for plotting

% Setup figure
figure(2)
set(gcf,'units','cent','position',[40,40,8,15])
clf

% Top panel: Ricker
subplot(2,1,1)
hold on
for x = 1:length(a)
Nt = a(x)*N.*exp(-b*N);
plot(N,Nt,'color','k','linestyle',LS{x})
end

% format plots
set(gca,'ylim',[0 18],'ytick',0:5:20,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:20:100)
ylabel(gca,'N_t_+_1','fontsize',12)
xlabel(gca,'N_t','fontsize',12)
text(2,17,'a','fontsize',12)

% Bottom panel: BH
subplot(2,1,2)
hold on
for x = 1:length(a)
Nt = N./(1/a(x) + N/b2);
plot(N,Nt,'color','k','linestyle',LS{x})
end

% format plots
set(gca,'ylim',[0 10],'ytick',0:2:20,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:20:100)
ylabel(gca,'N_t_+_1','fontsize',12)
xlabel(gca,'N_t','fontsize',12)
text(2,9.5,'b','fontsize',12)
