function totoaba

% Analysis of Totoaba model, based on Cisneros-Mata et al. (1995a, Ecol Appl)

A = 25; % 25 year classes
S = [0.01,0.01 repmat(0.764,1,22)]; % length A-1 survival vector 
Linf = 169.9; k = 0.152; t0 = -0.61; % von Bertalanffy parameters
Len = Linf*(1 - exp(-k*((1:A)-t0))); % length-at-age, from Cisneros-Mata et al. (1995b, Con Bio)
F = Len.^3; % assume fecundity proportional to length
F(1:6) = 0; % immature ages

% Assemble Leslie matrix:
L = [diag(S), zeros(length(S),1)];
L = [F; L];

% Rescale fecundity so that lambda = 1
l = max(eig(L)); % lambda
d = abs( l - 1); % difference between lambda & 1
%d = 0;
while d > 1e-4
    if l > 1
        L(1,:) = L(1,:)*0.99; % downscale
    else
        L(1,:) = L(1,:)*1.01; % downscale
    end
    l = max(eig(L)); % lambda
    d = abs( l - 1); % difference between lambda & 1
end

% Elasticities:
% first get the dominant eigenvectors:
[W,Ls,V]=eig(L);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
V = V(:,maxL);

% fecundity elasticities:
eF = V(1).*W'.*L(1,:)./l./dot(V,W);

% survival elasticities:
eS = V(2:end).*W(1:end-1).*S(:)./l./dot(V,W);

E = [eF(:);eS(:)];

% Project the model:
Lfec = L;
Lfec(1,:) = Lfec(1,:).*1.10; % 25% increase in fecundity
Lsurv = L;
Lsurv(2:6,:) = Lsurv(2:6,:)*1.10; % 25% increase in juvenile survival

T = 50;
Nf = nan(A,T);
N0 = W./W(end)*100; % initial population has 100 age-25 individuals
Nf(:,1) = N0;
Ns = Nf;
for t = 2:T
    Nf(:,t) = Lfec*Nf(:,t-1);
    Ns(:,t) = Lsurv*Ns(:,t-1);
end



% Figure:
% Top panel: survival & fecundity vs. age
figure(1)
clf
set(gcf,'units','cent','position',[2 30 9 21])

subplot(3,1,1)
set(gca,'position',[0.15 0.71 0.64 0.21])
cumS = [1,cumprod(S)];
hold on
pS = plot(1:A,log10(cumS),'k-');
set(gca,'color','none',...
        'ytick',[-8:2:0],'ylim',[-8 0],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[1 25],'xtick',[1 5 10 15 20])

ylabel('Survival to age (s_a)','fontsize',12)
axes('position',get(gca,'position'))
pA = plot(1:A,L(1,:),'k--');
set(gca,'yaxislocation','right',...
        'color','none',...
        'ytick',[0:2e3:10e3],'yticklabel',0:2:10,...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xlim',[1 25])
ylabel({'Fecundity at age (m_a)';'(thousands of age-1 females per female)'},'fontsize',12)
xlabel('Age (a)','fontsize',12)

subplot(3,1,2)
set(gca,'position',[0.15 0.41 0.64 0.21])
hold on
bar(E,'facecolor',[0.6 0.6 0.6])
set(gca,'color','none',...
        'ytick',0:0.02:0.1,'ylim',[0 0.1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xtick',[7 15 25 25+[10 20]],'xticklabel',[7 15 25 10 20],'xlim',[1 49])
ylabel('Elasticity','fontsize',12)
axes('position',get(gca,'position'))
plot(cumsum(E),'k-')
set(gca,'yaxislocation','right',...
        'color','none',... 
        'ytick',0:0.2:1,'ylim',[0 1],...
        'xtick',[],...
        'xlim',[1 49])
ylabel('Cumulative elasticity','fontsize',12)

subplot(3,1,3)
set(gca,'position',[0.15 0.11 0.64 0.21])
hold on
plot(1:T,sum(Nf(7:end,:)),'k--')
plot(1:T,sum(Ns(7:end,:)),'k-')
set(gca,'color','none',...
        'ytick',0:0.2e5:3e5,'ylim',[0 1e5],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xtick',0:5:30,'xlim',[1 30])
ylabel('Total population size','fontsize',12)
box on


% For Chapter 8, calculate the long-term variance in log(N) under
% different scenarios of environmental variability.
CV = 0:0.01:0.8;
Sig = nan(length(CV),6);
% block diagonal matrix for correlations among each set of parameters
Dm = [ones(25),zeros(25,24);...
      zeros(2,25),ones(2,2),zeros(2,22);...
      zeros(4,27),ones(4,4),zeros(4,18);...
      zeros(18,31),ones(18,18)];

  
for c = 1:length(CV)
% 1: r = 1, all parameters vary
D = ones(length(E))*CV(c)^2;
Sig(c,1) = E(:)'*D*E(:);

% 2 r = 1 only within parameter sets
Sig(c,2) = E(:)'*Dm*E(:)*CV(c)^2;

% 3: r = 0, only pre-adult survival
D = ones(4)*CV(c)^2;
Etmp = E((length(eF)+3):(length(eF)+6));
Sig(c,3) = Etmp(:)'*D*Etmp(:);

% 4: r = 0, only adult survival
D = ones(18)*CV(c)^2;
Etmp = E(32:end);
Sig(c,4) = Etmp(:)'*D*Etmp(:);

% 5: r = 0, only juvenile survival
D = ones(2)*CV(c)^2;
Etmp = E((length(eF)+1):(length(eF)+2));
Sig(c,5) = Etmp(:)'*D*Etmp(:);

% 6: r = 0, only reproduction
D = ones(25)*CV(c)^2;
Etmp = E(1:25);
Sig(c,6) = Etmp(:)'*D*Etmp(:);

% 7: r = 0, all parameters vary independently
D = eye(length(E))*CV(c)^2;
Sig(c,7) = E(:)'*D*E(:);
end

figure(2)
clf
set(gcf,'units','cent','position',[5 30 9 18])

% Plot elasticities again
subplot(2,1,1)
set(gca,'position',[0.13 0.58 0.7 0.34])
hold on
bar(E,'facecolor',[0.6 0.6 0.6])
set(gca,'color','none',...
        'ytick',0:0.02:0.1,'ylim',[0 0.1],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xtick',[7 15 25 25+[10 20]],'xticklabel',[7 15 25 10 20],'xlim',[1 49])
ylabel('Elasticity','fontsize',12)
axes('position',get(gca,'position'))
plot(cumsum(E),'k-')
set(gca,'yaxislocation','right',...
        'color','none',... 
        'ytick',0:0.2:1,'ylim',[0 1],...
        'xtick',[],...
        'xlim',[1 49])
ylabel('Cumulative elasticity','fontsize',12)

% Plot with full-scale y-axis
subplot(2,1,2)
hold on
Col=repmat(linspace(0,0.8,7)',[1,3]); % colormap
for c = [1 2 5]
%plot(CV.^2,Sig(:,c),'color',Col(c,:),'linewidth',1)
plot(CV,Sig(:,c),'color',Col(c,:),'linewidth',1)
end
set(gca,'ytick',0:0.05:1,'ylim',[0 0.25],...
        'tickdir','out','ticklength',[0.02 0.02],...
        'xtick',0:0.1:1,'xlim',[0 0.5],...
        'fontsize',10)
    Pos = get(gca,'position');
set(gca,'position',[Pos(1:2) 0.7 0.34])
ylabel('Rate of increase in variance (s2)','fontsize',12)
xlabel('CV2 of environmental variation in model parameters','fontsize',12)

