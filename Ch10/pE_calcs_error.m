function pE_calcs_error

% Create figure 10.1 showing the influence of estimation error on
% calculations of p(extinction). Based on Ludwig (1999) Table 1, Fieberg &
% Ellner (2000) Fig 1

% 1) Do the model simulations
n1 = 1e3; % number of sims for panel A
T = 1e2; % time horizon
T2 = 20; % shorter horizon
rs = -0.1:0.02:0.1; % possible true values of r
sigr = 0.2; % standard deviation in r
thresh = log(0.1); % extinction threshold (log of 10% of N0)
rr = normrnd(repmat(rs(:),[1,n1,T]),sigr);
sumr = cumsum(rr,3); % add up r's over time
E = mean(any(sumr<=thresh,3),2);

% 2) Diffusion approx
sigs = sigr;
mus = linspace(-0.1,0.1,1e2);
diffE = normcdf((thresh-mus.*T)./(sigs.*sqrt(T))) + exp(2.*mus.*thresh./(sigs.^2)).*normcdf((thresh+mus.*T)./(sigs*sqrt(T)));
diffE2 = normcdf((thresh-mus.*T2)./(sigs.*sqrt(T2))) + exp(2.*mus.*thresh./(sigs.^2)).*normcdf((thresh+mus.*T2)./(sigs*sqrt(T2)));


% 3) Simulate observations
Sampr = normrnd(0,0.2,[1e3,5]); % simulate samples of r with mu = 0, sig = 0.2
Meanr = mean(Sampr,2); % observed mean r
Stdr = std(Sampr,[],2); % observed sd of r
Simr = normrnd(repmat(Meanr,[1,n1,T]),repmat(Stdr,[1,n1,T]));
sumr = cumsum(Simr,3); % add up r's over time
Esim = mean(any(sumr<=thresh,3),2);

% 4) Simulation with shorter T, same true p[E]
Sampr2 = normrnd(-0.06,0.2,[100,5]); % simulate samples of r with mu = 0, sig = 0.2
Meanr2 = mean(Sampr2,2); % observed mean r
Stdr2 = std(Sampr2,[],2); % observed sd of r
Simr2 = normrnd(repmat(Meanr2,[1,n1,T2]),repmat(Stdr2,[1,n1,T2]));
sumr2 = cumsum(Simr2,3); % add up r's over time
Esim2 = mean(any(sumr2<=thresh,3),2);


figure(1)
clf
set(gcf,'units','cent','position',[10,10,21,24])

% Plot 1: p(E) vs r
subplot(2,2,1)
hold on
plot(rs,E,'ko')
plot(mus,diffE,'k')
%plot(mus,diffE2,'r')
plot([-0.2 0.2],[0.2 0.2],'k--')
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',12)
xlabel('Intrinsic growth rate, r','fontsize',14)
ylabel('Probability of extinction, p(E)','fontsize',14)

% Plot 2: observed mean r's
subplot(2,2,3)
hold on
plot(Meanr,Esim,'ko');
set(gca,'xlim',[-0.2 0.2])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',12)
xlabel('Estimated value of mean r','fontsize',14)
ylabel('Probability of extinction, p(E)','fontsize',14)

% inset - histogram of r
axes('position',[0.35 0.32 0.1 0.1])
histogram(Meanr,10,'norm','prob','facecolor',[0.5 0.5 0.5])
set(gca,'ytick',[])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'xlim',[-0.2 0.2])


% Plot 3: p(E)
subplot(2,2,2)
hold on
%histogram(Esim2,10,'norm','prob','facecolor','r')
histogram(Esim,10,'norm','prob','facecolor',[0.5 0.5 0.5])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',12)
xlabel('Probability of extinction, p(E)','fontsize',14)
ylabel('Frequency','fontsize',14)


% Plot 4: data from Ludwig (1999)
Phi = [0 0.25 0.5];
Lud_r(:,:,1) = [0.33, 1.44; 0.089, 1.65; -0.75, 2.36];
Lud_e(:,:,1) = [0, 0.035; 0 0.4; 0, 1];

Lud_r(:,:,2) = [-0.082, 0.64; -0.12, 0.5; -0.11, 0.4];
Lud_e(:,:,2) = [0, 0.99; 0 1; 0, 1];

Lud_r(:,:,3) = [-0.66, 1.7; -0.92, 1.74; -1.02, 1.63];
Lud_e(:,:,3) = [7.4e-6, 0.035; 0 0.4; 0, 1];

Pos = [0.11, 0.22, 0.33];

Lims = [-5.4, 5.4; -1.2 1.2; -2.3, 2.3];

for i = 1:3

axes('position',[0.5703, Pos(i), 0.15 0.1])
hold on
plot(Phi,Lud_r(:,1,i),'ko')
plot(Phi,Lud_r(:,2,i),'kd')
set(gca,'xlim',[-0.1 0.6])
if i > 1; set(gca,'xtick',[]); end
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[Lims(i,1),Lims(i,2)])

axes('position',[0.73, Pos(i), 0.15 0.1],'YAxisLocation','right')
hold on
plot(Phi,Lud_e(:,1,i),'k')
plot(Phi,Lud_e(:,2,i),'k')
set(gca,'ylim',[-0.1 1.1],'color','none','xlim',[-0.1 0.6])
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])

if i > 1; set(gca,'xtick',[]); end

end


function prettify
set(gca,'xcolor','k','ycolor','k','tickdir','out','ticklength',[0.015 0.015])


