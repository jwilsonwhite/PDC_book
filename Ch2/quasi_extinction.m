function quasi_extinction

% Plot an example of a discrete-time stochastic model that is close to
% extinction

Lm = 1; % mean lambda
Ls = 0.2; % sd lambda
X = 1e4; % Number of simulated populations
T = 30; % time horizon for extinction
Nmin = 4; % quasi-extinction threshold (40% of initial abundance)

% Pre-allocate state variable
N = ones(T,X)*10;

% Create matrix of growth rates
L = normrnd(Lm,Ls,[T,X]);

% Iterate model:
for t = 2:T
    N(t,:) = N(t-1,:).*L(t,:);
end

% Calculate & print some summary statistics
mean(N(end,:))
exp(mean(log(N(end,:))))
10*Lm^T
exp((log(Lm) - (Ls^2)/(2*Lm^2))*T)*10


% Get histogram of final abundances:
Edges = [0:1:30,Inf]; % bins
H = histc(N(end,:),Edges); % histogram
H(end-1) = H(end-1)+H(end); % pool last bin

% Plot a few example trajectories, plus histogram
figure(1)
clf
set(gcf,'units','cent','position',[10,40,16,8])

subplot(1,3,1:2)
hold on
set(gca,'position',[0.13, 0.11, 0.56,0.815])
Pos = N(end,:)>Nmin;
NPos = N(:,Pos);
Neg = N(end,:)<Nmin;
NNeg = N(:,Neg);
ph = plot(0:T-1,NPos(:,1:4),'k');
ph = plot(0:T-1,NNeg(:,1),'color',[0.5 0.5 0.5]);
plot([0,T-1],[Nmin, Nmin],'k--')
%IsExt = N(end,1:5) < Nmin;
%set(ph(IsExt),'color',[0.5 0.5 0.5])
ylim([0 30])
xlim([0 T-1])
set(gca,'xtick',0:5:25)
xlabel('Time','fontsize',12)
ylabel('Population size (N_t)','fontsize',12)


subplot(1,3,3)
set(gca,'position',[0.69, 0.11, 0.24,0.815])
hold on
barh(Edges(1:Nmin)+0.5,H(1:Nmin),1,'facecolor',[0.5 0.5 0.5])
barh(Edges((Nmin+1):end-1)+0.5,H((Nmin+1):end-1),1,'facecolor',[1 1 1])
shading flat
ylim([0 30])
set(gca,'ytick',[],'xtick',[])
text(-100,-3,'Proportion of simulations','fontsize',12)
%N2(N2>5) = NaN;
%hist(N2,100)