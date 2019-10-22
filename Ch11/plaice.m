function plaice

% North Sea Plaice model of Beverton and Holt (1957)

A = 15; % max age
Ages = 2:A;
A0 = min(Ages);
R  = 1; % mean recruitment
M = 0.1; % natural mortality
Winf = 2867; % asymptotic maximum weight
K = 0.095; % von Bert growth rate
t0 = -0.815; % time at L = 0
Linf = 68.5; % von Bert asymptotic maximum length
L = Linf .* (1 - exp(-K.*(Ages-t0))); % vB length
W = Winf .* (1 - exp(-K.*([Ages,A+1]-t0))).^3; % vB weight
W = W(2:end);

As = 4:1:15; % range of ages of entry to fishery
Fs = [0:0.01:1.5, 10]; % range of possible Fs

% vectorize to create grids for evaluation
Asmat = repmat(As(:),[1 length(Fs)]);
Asmat = repmat(Asmat,[1,1,length(Ages)]);
Fsmat = repmat(Fs(:)',[length(As),1]);
Fsmat = repmat(Fsmat,[1,1,length(Ages)]);
Agemat = repmat( reshape(Ages,[1,1,length(Ages)]),[length(As),length(Fs),1]);
Wmat = repmat( reshape(W,[1,1,length(Ages)]),[length(As),length(Fs),1]);

% Calculate equilibrium population age distribution
N = R.* exp( -M.*(Agemat-A0)) .*exp( -Fsmat.*(max(0,Agemat-Asmat)));
% Number of dead in each age class
Dead = N(:,:,1:end-1) - N(:,:,2:end);
Dead(:,:,end+1)=0 ;%N(:,:,end);
% Proportion of dead fish that were harvested, converted to biomass yield
Yield = Dead.*Fsmat./(Fsmat+M).*...
    (Agemat>=Asmat)...
    .*Wmat;

Y = sum(Yield,3);

% Max yield, find F_msy for each age of entry
MaxY = max(Y,[],2);
MaxAs = Y==repmat(MaxY(:),[1,size(Y,2)]);
MaxAs(end,:) = false; % take off last row, all zeros
Asmat2 = Asmat(:,:,1);
Fsmat2 = Fsmat(:,:,1);

MA = Asmat2((MaxAs));
MF = Fsmat2((MaxAs));


% Plotting the contours:
figure(1)
clf
set(gcf,'units','cent','position',[5 30 9 15])

subplot(2,1,1)
hold on

colormap(1-repmat((0.2:0.05:0.7)',[1,3]))
Fplot = Fs;
Fplot(end) = Fplot(end-1) + 10*diff(Fplot(1:2));
contourf(Fplot,As,Y,0:50:500)

set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',0:0.5:1.5,'ytick',4:2:16)
set(gca,'xcolor','k','ycolor','k')
axis square

ch = colorbar;
set(ch,'tickdir','out','ticklength',0.015)
set(ch,'ycolor','k')

MF(MF>1.5) = Fplot(end);
plot(MF,MA,'ko','color','k','linestyle','-')

subplot(2,1,2)
hold on

contourf(Fplot,As,max(0,Y - 100.*Fsmat2),0:50:500)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',0:0.5:1.5,'ytick',4:2:16)
set(gca,'xcolor','k','ycolor','k')
axis square

ch = colorbar;
set(ch,'tickdir','out','ticklength',0.015)
set(ch,'ycolor','k')

%plot(MF,MA,'ko','color','k','linestyle','-')
