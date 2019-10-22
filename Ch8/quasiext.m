function quasiext

% quasi-extinction probabilities - example from Cisneros-Mata et al. (1997)
% using the Lande & Orzack formula

% variances (sigma^2)
Sig2 = 0.001:0.001:0.5;
Mu = -Sig2/2;
Ts = [5 10 25 50 75 100]; % time horizons
Xs = [0.05 0.1 0.2 0.3 0.4 0.5]; % extinction thresholds (as % reductions)

for t = 1:length(Ts)
    for x = 1:length(Xs)
       
        % this is the way that Lande & Orzack wrote it, rearranged from the
        % way Cisneros-Mata & Botsford did
     
        Px(:,x,t) = normcdf( (-log(1/Xs(x)) - Mu.*Ts(t))./sqrt(Sig2.*Ts(t)))...
             + exp((-2.*Mu.*log(1/Xs(x)))./Sig2).*...
             (1-normcdf( (log(1/Xs(x)) - Mu.*Ts(t))./sqrt(Sig2.*Ts(t))));
        
    end
end

figure(3)
clf
set(gcf,'units','cent','position',[10 10 9 16])



Col = 1-gray(10);
Col = Col(4:end,:);
Lty = {'-.','--','-'};

% Panel 1: time horizon varies
sh(1)=subplot(2,1,1);
hold on
X = 2; % hold threshold constant
for t = 1:length(Ts)
plot(Sig2,Px(:,X,t),'color',Col(t,:),'linestyle','-');
end
  
% Panel 2: threshold varies
sh(2)=subplot(2,1,2);
hold on
T = 3; % hold threshold constant
for x = 1:length(Xs)
plot(Sig2,Px(:,x,T),'color',Col(x,:),'linestyle','-');
end

set(sh(:),'tickdir','out','ticklength',[0.015 0.015],...
         'xcolor',zeros(3,1),'ycolor',zeros(3,1),...
         'ylim',[0 1],'ytick',0:0.2:1,'xlim',[0 0.5],'xtick',0:0.1:0.5,'fontsize',10)
 ylabel(sh(1),'Probability of quasi-extinction by time t, pt(X|N0)','fontsize',12)
 xlabel(sh(2),'Rate of increase in variance (s2)','fontsize',12)

