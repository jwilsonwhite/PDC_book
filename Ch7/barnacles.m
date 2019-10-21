function barnacles

% Implement Roughgarden et al. (1985, Ecology) model of barnacle dynamics

% Parameters from the paper:
Ages =   1:103; % 103 weeks
A = 100; % 100 cm2 total area
P0 = 0.9; % per week
m = 35; % density-dependence control parameter
Cexp = 4; % growth exponent
C0 = 7.07/(103^Cexp); % growth coefficient, based on 7.07cm2 at 103 weeks
Ca = C0.*Ages.^Cexp; % size at age

T = 1000; % 1000 weeks
S = [5, 15, 25]; % settlement rates (one row of figure for each rate. Higher settlement = lower stability)

figure(1)
clf
set(gcf,'units','cent','position',[10 10 15 10])

Sp1 = [1 4 7];
Sp2 = [2 5 8];
Sp3 = [3 6 9];

for s = 1:length(S)
    
    % Initialize variables
    F = zeros(T,1);
    F(1) = A;
    N = zeros(length(Ages),T);
    
    % Loop over time
    for t = 2:T
        
        R = S(s)*F(t-1); % settlement (depends on availability of free space)
        N(2:end,t) = P0*N(1:end-1,t-1); % adult survival
        N(1,t) = R;
        
        F(t) = A - Ca(:)'*N(:,t); % recalculate free space
        F(t) = max(0,F(t));
  
    end
    
    % alternative version including density-dependent mortality, too
    F2 = zeros(T,1);
    F2(1) = A;
    N = zeros(length(Ages),T);
    for t = 2:T
        
        R = S(s)*F(t-1);
        Px = P0*(1-exp(-m*F2(t-1)/A));
        N(2:end,t) = Px*N(1:end-1,t-1);
        N(1,t) = R;
        
        F2(t) = A - Ca(:)'*N(:,t);
        F2(t) = max(0,F2(t));
  
    end
    
       % Calculate K (the elasticity of the density-dependent function)...
       Surv = P0.^(Ages-1);
    PhiC = Surv.*Ca;
Rstar = S(s)*A/(1+S(s)*sum(PhiC));
K = -1*Rstar*sum(PhiC)/(A-sum(PhiC)*Rstar);
disp(K)
    
    % Plot influence function
    subplot(3,3,Sp1(s))
    hold on
    
    plot(Ages,K*PhiC./sum(PhiC),'k-');
    plot(Ages,zeros(length(Ages),1),'k-')
    plot([Ages(end),Ages(end)],[0 K*PhiC(end)/sum(PhiC)],'k-')
    
    set(gca,'tickdir','out','ticklength',[0.015 0.015])
    set(gca,'xtick',1:25:100,'xticklabel',0:0.5:2,'ytick',0:0.2:1)
    set(gca,'fontsize',10)
    set(gca,'ylim',[-0.08 0.03],'xlim',[0 115])
    
    % Plot f(S) & equilbrium line
    subplot(3,3,Sp2(s))
    hold on
    
    Rs = 0:1e3;
   plot(Rs,S(s).*(A-sum(PhiC).*Rs(:)),'k-')
   plot([0,max(Rs)],[0 max(Rs)/sum(PhiC)],'k--')
    set(gca,'tickdir','out','ticklength',[0.015 0.015])
    set(gca,'xtick',[],'ytick',[])
    set(gca,'fontsize',10)
    set(gca,'xlim',[0 max(Rs)*0.9],'ylim',[0 S(s)*A*1.05])
    
    
    % Plot timeseries
    subplot(3,3,Sp3(s))
    hold on
    plot(1:T,F/A,'k-')
    plot(1:T,F2/A,'k-','color',[0.5 0.5 0.5])
    
    set(gca,'tickdir','out','ticklength',[0.015 0.015])
    set(gca,'xtick',1:104:2000,'xticklabel',0:2:520,'ytick',0:0.5:1)
    set(gca,'fontsize',10)
    set(gca,'ylim',[0 1])
    
end
