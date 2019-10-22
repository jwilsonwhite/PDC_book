function fox

% Figure 10.5 - Bakker et al's island fox model

% Panel 1: Survival as a function of adult density:

% Regression parameters from Table 3:
b = [0.586, -0.464, -0.523, -12.041, 60.7, -5, -0.72 1.215, -0.931]; % NB it appears the digits were transposed for pup effect in printed table (assuming their graph is correct)
bsd = [0.284, 0.12, 0.162, 20.06, 11.8, 0.982, 0.327, 0.42, 0.548];

% Breeding params (Table 6):
B1 = 0.429; B1sd = 0.089; % pups (Santa Cruz)
B2 = 0.618; B2sd = 0.124; % adults (Santa Cruz)

% Litter size (Table 6):
L = 1.92; Lsd = 0.119;


Dens = linspace(0,15,100);

% Two initial examples: adults, pups, + 13 eagles
X = b(1) + b(5)*(Dens/100) + b(6)*(Dens.^2)/100;
Surv1 = 1./(1+exp(-X));

Xpup = b(1) + b(3) + b(5)*(Dens/100) + b(6)*(Dens.^2)/100;
Survpup1 = 1./(1+exp(-Xpup));

XE = b(1) + b(5)*(Dens/100) + b(6)*(Dens.^2)/100 + b(4)*0.13;
SurvE = 1./(1+exp(-XE));

XpupE = b(1) + b(3) + b(5)*(Dens/100) + b(6)*(Dens.^2)/100 + b(4)*0.13;
SurvpupE = 1./(1+exp(-XpupE));



%------------------------------
% Plotting
figure(1)
clf
set(gcf,'units','cent','position',[10 10 20 20])

% Panel 1: survival
subplot(3,2,1)
hold on
plot(Dens,Surv1,'k','linewidth',1)
plot(Dens,Survpup1,'color',[0.5 0.5 0.5],'linewidth',1)
plot(Dens,SurvE,'k--','linewidth',1)
plot(Dens,SurvpupE,'k--','color',[0.5 0.5 0.5],'linewidth',1)
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
set(gca,'ytick',0:0.2:1)
xlabel('Adult fox density (km-2)','fontsize',12)
ylabel('Annual survival','fontsize',12)
%------------------------------

% Get deterministic projection matrix:
for i = 1:length(Surv1)

    M = [0.5*Survpup1(i)*B1*L, 0.5*Surv1(i)*B2*L, 0, 0;...
        Survpup1(i), Surv1(i), 0, 0;...
        0.5*Survpup1(i)*B1*L, 0.5*Surv1(i)*B2*L, 0, 0;...
        0, 0, Survpup1(i), Surv1(i)];
    Eig(i) = max(eig(M));

end

%-------------------------
% Stock-recruit plot
subplot(3,2,2)
hold on
plot(Dens,Dens.*Surv1.*Eig,'k')
plot(Dens,Dens,'k--')
set(gca,'tickdir','out','ticklength',[0.02 0.02])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
xlabel('Adult fox density, time t (km-2)','fontsize',12)
ylabel('Adult fox density, time t+1 (km-2)','fontsize',12)
%--------------------------



% Make some simulated runs:
T = 50;
NN = 10; % number of sims
Dens_conv = 250; % this converts density to total numbers (total area of Santa Cruz)

n0 = 650*[0.35/2 0.35/2 0.65/2 0.65/2]'; % initial density + stage distribution
n0 = round(n0);

% ENSO has 11% prob every year (binary)
ENSO = unifrnd(0,1,NN,T)<= 0.11;

% How many eagles
Eagle = 6;

for nn = 1:NN
    
    N1(:,1,nn) = n0; % baseline
    NE(:,1,nn) = n0; % plus eagles
    NU(:,1,nn) = n0; % uncertainty
    NUE(:,1,nn) = n0; % uncertainty plus eagles
    
    for t = 2:T
        
    % Baseline: 
    N1(:,t,nn) = run_model(N1(:,t-1,nn),Dens_conv,b,0,ENSO(nn,t-1:t),0,B1,0,B2,0,L,0);

    % Eagles:
    NE(:,t,nn) = run_model(NE(:,t-1,nn),Dens_conv,b,0,ENSO(nn,t-1:t),Eagle,B1,0,B2,0,L,0);
    
    % Baseline: 
    NU(:,t,nn) = run_model(NU(:,t-1,nn),Dens_conv,b,bsd,ENSO(nn,t-1:t),0,B1,B1sd,B2,B2sd,L,Lsd);

    % Eagles:
    NUE(:,t,nn) = run_model(NUE(:,t-1,nn),Dens_conv,b,bsd,ENSO(nn,t-1:t),Eagle,B1,B1sd,B2,B2sd,L,Lsd);

    
        
    end % end T
end % end NN


%-------------------------
% Plot time series
subplot(3,2,3)
hold on
N = squeeze(sum(N1(3:4,:,:),1));
isExt = N(end,:)<=30;
plot(N(:,isExt),'k')
plot(N(:,~isExt),'color',repmat(0.7,[1,3]))

set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
set(gca,'ylim',[-200 6000])
xlabel('Time (y)','fontsize',12)
ylabel('Adult fox abundance','fontsize',12)
%--------------------------

%-------------------------
% Plot time series
subplot(3,2,4)
hold on
N = squeeze(sum(NE(3:4,:,:),1));
isExt = N(end,:)<=30;
plot(N(:,isExt),'k')
plot(N(:,~isExt),'color',repmat(0.7,[1,3]))

set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
set(gca,'ylim',[-200 6000])
xlabel('Time (y)','fontsize',12)
ylabel('Adult fox abundance','fontsize',12)
%--------------------------


%-------------------------
% Plot time series
subplot(3,2,5)
hold on
N = squeeze(sum(NU(3:4,:,:),1));
isExt = N(end,:)<=30;
plot(N(:,isExt),'k')
plot(N(:,~isExt),'color',repmat(0.7,[1,3]))

set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
set(gca,'ylim',[-200 6000])
xlabel('Time (y)','fontsize',12)
ylabel('Adult fox abundance','fontsize',12)
%--------------------------


%-------------------------
% Plot time series
subplot(3,2,6)
hold on
N = squeeze(sum(NUE(3:4,:,:),1));
isExt = N(end,:)<=30;
plot(N(:,isExt),'k')
plot(N(:,~isExt),'color',repmat(0.7,[1,3]))

set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'xcolor','k','ycolor','k','fontsize',10)
set(gca,'ylim',[-200 6000])
xlabel('Time (y)','fontsize',12)
ylabel('Adult fox abundance','fontsize',12)
%--------------------------


function S= get_surv(b,bse,Dens,Eagle,ENSO,Pup)
% calculate survival rate
b = normrnd(b,bse);

X = b(1) + b(3)*Pup + b(4)*Eagle/100 + b(5)*Dens/100 + b(6)*(Dens^2)/100 + b(7)*ENSO(1) + b(8)*ENSO(2)*Pup;
Surv = 1/(1+exp(-X));

if Pup
    Vu = 0.094;
else
Vu = 0.046; % overall variance in survival
end
Var = Vu*Surv*(1-Surv); % variance expression used by Bakker
% params of beta:
aa = ((1 - Surv)/Var - 1/Surv)*Surv^2;
bb = aa*(1/Surv - 1);
S = betarnd(aa,bb);


function B = get_birth(B1,B1se,B2,B2se,L,Lse,N)

B1t = min(1,max(0,normrnd(B1,B1se)));
B2t = min(1,max(0,normrnd(B2,B2se)));
L = max(0,normrnd(L,Lse));

% Adjust for imbalanced sex ratios
B2tt = min(B2t,B2t*N(3)/N(4)); % penalize if too few males
B1tt = min(1,B1t*(N(3) - B2tt*N(4))/(N(2))); % adult females get the males first

% total births
%Btotal = round( L*binornd(N(4),B2tt) + L*binornd(N(2),B1tt));
Btotal = poissrnd(L*N(4)*B2tt) + poissrnd(L*N(2)*B1tt);

% randomly assign to sex
B(1) = binornd(Btotal,0.5);
B(2) = Btotal-B(1);

if any(isnan(B)); keyboard; end

function N_out = run_model(N_in,Dens_conv,b,bse,ENSO,Eagle,B1,B1se,B2,B2se,L,Lse)

    Dens = sum(N_in(3:4))/Dens_conv;
    Surv = get_surv(b,bse,Dens,Eagle,ENSO,0);
    Survpup = get_surv(b,bse,Dens,Eagle,ENSO,1);
    
    % Advance survival, apply demographic stochasticity
    N_out(3:4) = binornd(N_in(3:4),Surv) + binornd(N_in(1:2),Survpup);  
    
    % Reproduction
    N_out(1:2) = get_birth(B1,B1se,B2,B2se,L,Lse,N_in(:));
    
    % Quasi-extinction
    if sum(N_out(3:4))<=30; N_out= 0*N_out; end









