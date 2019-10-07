function theta_logistic_plots

% Create example plots for theta logistic model

% Same parameters used in the logistic plots earlier
R = [exp(0.4) exp(0.2) exp(0.1)]; % growth rates
K =10; % carrying capacity
Th = [1/3 1 3]; % theta values
LS = {':','-','--'}; % line types for plotting


% Setup figure
figure(2)
set(gcf,'units','cent','position',[40,40,8,15])
clf

% Upper panel: dN/dt * 1/N
subplot(2,1,1)
hold on
N = 0:0.01:K; % values of N
for t = 1:3
% get theta logistic solution
Nt = Theta_logistic(1,N,log(R(t)),K,Th(t));
% plot results
ph = plot(N,Nt./N,'color','k','linestyle',LS{t});
end

% Format plots
set(gca,'ytick',0:0.1:2,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:5:10,'xticklabel',{'0','K/2','K'})
ylabel(gca,'Per capita growth rate (dN/dt 1/N)','fontsize',12)
xlabel(gca,'Population size (N)','fontsize',12)
text(0.23,0.38,'a','fontsize',12)

% Lower panel: dN/dt
subplot(2,1,2)
hold on
N = 0:0.01:K;
for t = 1:3
    % get theta logistic solution
Nt = Theta_logistic(1,N,log(R(t)),K,Th(t));
% plot results
ph = plot(N,Nt,'color','k','linestyle',LS{t});
end

% Format plots
set(gca,'ytick',0:0.1:2,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:5:10,'xticklabel',{'0','K/2','K'})
ylabel(gca,'Overall growth rate (dN/dt)','fontsize',12)
xlabel(gca,'Population size (N)','fontsize',12)
text(0.23,0.48,'b','fontsize',12)


% Define model functions
function dN = Theta_logistic(t,N,r,K,Th)
% Exponential growth
dN = r.*N.*(1-(N./K).^Th);