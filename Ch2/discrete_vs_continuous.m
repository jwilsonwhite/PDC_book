function discrete_vs_continuous

% Plot examples of exponential (continuous time) & geometric (discrete
% time) growth & decay

R = [1.15, 0.8]; % positive and negative growth rates

% time span over which to solve model
T = 0:10;
tspan = [T(1), T(end)];

% initial conditions
N0 = 1;

% Setup figure
figure(1)
set(gcf,'units','cent','position',[40,40,8,15])
clf

% Continuous:
sh(1) = subplot(2,1,1);
hold on

% Solve the ODE for each value of R:
for r = 1:length(R)
[TT,N] = ode23(@(t,N) ExpGrow(t,N,log(R(r))),tspan,N0);
ph=plot(TT,N,'k-');
end

% Discrete:
sh(2) = subplot(2,1,2);
hold on

% Iterate the discrete model for each value of R
for r = 1:length(R)
    Nt = N0.*R(r).^T;
    ph = plot(T,Nt,'ko');
end

% Plot cleanup:
set(sh(:),'tickdir','out','ticklength',[0.02 0.02])
xlabel(sh(2),'Time (y)','fontsize',12)
ylabel(sh(1),'Population size (N(t))','fontsize',12)
ylabel(sh(2),'Population size (N_t)','fontsize',12)
set(sh(:),'ytick',0:1:10,'xtick',0:2:10,'ylim',[0 4])

axes(sh(1))
text(0.2,3.7,'a','fontsize',12)
axes(sh(2))
text(0.2,3.7,'b','fontsize',12)

function dN = ExpGrow(t,N,r)
% Exponential growth
dN = r*N;