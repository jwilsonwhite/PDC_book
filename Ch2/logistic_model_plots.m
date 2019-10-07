function logistic_model_plots

% Create example plots for logistic model

R = 2; % growth rate
K = 10; % carrying capacity
N0 = 0.1; % initial conditions

tspan = [0 20]; % time interval for model solution

% Setup figure
figure(2)
set(gcf,'units','cent','position',[40,40,8,15])
clf

% Right panel: solution to model
subplot(1,2,2)
[TT,N] = ode45(@(t,N) Logistic(t,N,log(R),K),tspan,N0);
ph=plot(TT,N,'k-');
set(gca,'ylim',[0 11],'ytick',0:2:10,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:5:20)
ylabel(gca,'Population size (N)','fontsize',12)
xlabel(gca,'Time','fontsize',12)
text(0.5,10.2,'b','fontsize',12)

% Left panel: dN/dt
subplot(1,2,1)
N = 0:0.01:K;
Nt = Logistic(1,N,log(R),K);
ph = plot(N,Nt,'k-');
set(gca,'ylim',[0 2],'ytick',0:0.4:2,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xtick',0:2.5:10)
ylabel(gca,'Growth rate (dN/dt)','fontsize',12)
xlabel(gca,'Population size (N)','fontsize',12)
text(0.25,1.9,'a','fontsize',12)

% Define the model equation
function dN = Logistic(t,N,r,K)
% Exponential growth
dN = r.*N.*(1-N./K);