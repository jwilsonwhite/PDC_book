function chinook_extinct2

% Create Fig. 10.3b - example of PVA for Sacramento winter chinook with
% measurement error

% Model parameters:
% spawning age distribution
Sig2 = 0.25;
Sig3 = 0.67;
Sig4 = 0.08;
Sig = [Sig2 Sig3 Sig4];

SD_e = 1; % variance in cohort replacement rate
QE = 100; % quasi-extinction threshold

%Ninit = [1e3 2.5e3 5e3 7.5e3 1e4 2.5e4 5e4 7.5e4 1e5]; % (initial
%parameter exploration)
Ninit = [1e3:2e3:9e3 1e4:2e4:9e4 1e5]; % Initial popluation abundance
SigM = 0:0.25:1.5; % measurement error

E = 0;

N = 3:20; % number of samples 

T = 50; % time horizon
n = 5e4; % number of random simulations

% Loop over Ninit & Number of samples for each value of SigM to generate
% contour plots
for s = 1:length(SigM)
for ni = 1:length(Ninit)
for nn = 1:length(N)

% Mean replacement rate, as estimated from N(nn) samples with given measurement
% error
Emean = normrnd(E,sqrt(SD_e^2+2*SigM(s)^2)/sqrt(N(nn)-1),[n,1]);

% Generate distribution of replacement rates
Emat = exp(normrnd(repmat(Emean,[1,T]),SD_e,[n,T]));
    
% Initial spawner abundance
S1 = Ninit(ni);

% Pre-allocate state variables for simulation
S = zeros(n,T);
S(:,1:4) = S1;
    
    for t = 5:T
        
      % Check for spawning < QE
      Spawners = [S(:,t-2), S(:,t-3), S(:,t-4)];
      Repro = [Emat(:,t-2), Emat(:,t-3), Emat(:,t-4)];
      
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

   Ext(nn,ni,s) = mean(isExt);
    
end % end E loop
end % end Ninit loop
end % end SigM loop

% Plotting
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 8])
hold on

for s = 1:length(SigM)
contour(Ninit,N,Ext(:,:,s),[0.1,0.1],'color','k');
end
set(gca,'Xscale','log')
set(gca,'tickdir','out','ticklength',[0.03 0.03])
set(gca,'xcolor','k','ycolor','k')
axis square
