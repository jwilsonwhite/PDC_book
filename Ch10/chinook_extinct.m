function chinook_extinct

% Create Fig. 10.3 - example of PVA for Sacramento winter chinook

% Model parameters:
% spawning age distribution
Sig2 = 0.25;
Sig3 = 0.67;
Sig4 = 0.08;
Sig = [Sig2 Sig3 Sig4];

SD_e = 1; % variance in cohort replacement rate

QE = 100*2; % quasi-extinction threshold

Ninit = [1e3 2.5e3 5e3 7.5e3 1e4 2.5e4 5e4 7.5e4 1e5]; % initial population size
E = -0.6:0.05:0.1; % mean cohort replacement rate 


T = 50; % time horizon
n = 1e4; % number o simulations

for ni = 1:length(Ninit)
for e = 1:length(E)

    % random draws of cohort replacement rate
Emat = exp(normrnd(E(e),SD_e,[n,T]));
    
% Initial spawner abundance
S1 = Ninit(ni);

S = zeros(n,T);
S(:,1:4) = S1;
    
    for t = 5:T
        
    % Assemble variables to make model advancing simpler
      Spawners = [S(:,t-2), S(:,t-3), S(:,t-4)];
      Repro = [Emat(:,t-2), Emat(:,t-3), Emat(:,t-4)];

      % update age distribution
      S(:,t) = Sig2*Repro(:,1).*Spawners(:,1) + Sig3*Repro(:,2).*Spawners(:,2) + Sig4*Repro(:,3).*Spawners(:,3);
      
      % zero out repro if that spawning class is too small
      Emat(S(:,t)<QE,t) = 0;
        
    end % end T
    

% Count up extinctions over 3 temporal cohorts (bc extinction requires 3
% consecutive low years)
for j = 3:T

   S_ext = sum(S(:,j-2:j)<QE,2);

end


   isExt = any(S_ext == 3,2);

   Ext(e,ni) = mean(isExt);
    

end % end E loop
end % end Ninit loop

% Plotting
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 8])
hold on
contourf(Ninit,E,Ext,[0.1:0.1:1],'color','k');
contour(Ninit,E,Ext,[0.01 0.05],'color',[0 0 0]);
set(gca,'Xscale','log')
set(gca,'tickdir','out','ticklength',[0.03 0.03])
set(gca,'xcolor','k','ycolor','k')
axis square

G = 1-gray;
colormap(G,10)
ch = colorbar;
set(ch,'tickdir','out','ticklength',0.01)
