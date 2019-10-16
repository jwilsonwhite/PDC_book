function stage_projection

% Project a size-structured model and plot in terms of size. Designed to
% illustrate silliness of stage-structured models.

% Assemble the matrix (for loggerhead turtles, from Crowder et al. 1994):
F = [0 0 0 4.665 61.896]; % fecundity
D = [0 0.703 0.657 0.682 0.8091]; % diagonal (stage retention)
S = [0.675 0.047 0.019 0.061]; % survival
Smat = [diag(S), zeros(length(S),1)]; % lower half of matrix
L = [F;Smat]+diag(D);  % the full matrix

% Size-age relationship (from Crouse, citing Frazer):
Len = [5 mean([10.1 58]), mean([58.1 80]), mean([80.1 87]), 88];

% SAD:
N0 = zeros(size(L,1),1);
N0(1) = 1000;

% zero out reproduction:
Lmaster = L;
L(1,:) = 0;

% Project:
T = 54; % max age of the turtles
N = nan(length(N0),T);
N(:,1) = N0;
for t = 2:T
    N(:,t) = L*N(:,t-1);
end

SD = N;
% scale to create size distribution:
SD = SD./repmat(sum(SD),[size(L,1),1]);




%--------------------
% Elasticities:
% first get the dominant eigenvectors:
[W,Ls,V]=eig(Lmaster);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
V = V(:,maxL);
l = max(Ls);

% Stage durations:
Ts = [1 7 8 6 32];
% Survival
sig = [0.6747, 0.75, 0.6758, 0.7425, 0.8091];

% fecundity elasticities:
eF = V(1).*W'.*L(1,:)./l./dot(V,W);


% fecundity sensitivity:
sF = V(1).*W'./dot(V,W);

% growth elasticities:
eS = V(2:end).*W(1:end-1).*S(:)./l./dot(V,W);

% survival-within-stage sensitivity
sP = V(:).*W(:)./dot(V,W);

% growth senstivity
sG = V(2:end).*W(1:end-1)./dot(V,W);
sG(5) = 0;

% sensitivity of growth to survival:
sGsig = (sig./l).^(Ts-1)./(( (sig./l).^Ts - 1).^2) .* (Ts.*(1-sig./l) +  (sig./l).^(Ts+1) -(sig./l));
sGsig(5) = 0;
% sensitivity of within-stage-survival to survival
sPsig = 1 - sGsig;

% elasticity to sigma:
E = sig(:)./l.*(sP(:).*sPsig(:) + sG(:).*sGsig(:));


figure(2)
set(gcf,'units','cent','position',[5 5 9 18])
clf
subplot(2,1,1)
bar(E,0.8,'facecolor',[0.5 0.5 0.5])
set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'ylim',[0 0.5],'ytick',0:0.1:0.5)
set(gca,'fontsize',10)
xlabel('Stage','fontsize',12)
ylabel('Sensitivity of l to survival','fontsize',12)

%------------------------
% Transient model
Fec = [zeros(1,sum(Ts(1:3))), zeros(1,Ts(4))*F(4), ones(1,Ts(5))*76.5];
Surv = [ones(1,Ts(1))*sig(1),...
    ones(1,Ts(2))*sig(2),...
    ones(1,Ts(3))*sig(3),...
    ones(1,Ts(4))*sig(4),...
    ones(1,Ts(5))*sig(5)];
Surv = Surv(1:end-1);
LA = diag(Surv);
LA = [LA, zeros(length(Surv),1)];
LA = [Fec; LA];

% Do cohort advance:
T = 54; % max age of the turtles
N = nan(T,T);
N(:,1) = [1; zeros(T-1,1)] ;
LAt = LA;
LAt(1,:)= 0; % zero out reproduction
for t = 2:T
    N(:,t) = LAt*N(:,t-1);
end

SDa = zeros(5,T);
% translate back to stages:
SDa(1,:) = (N(1,:));
SDa(2,:) = sum(N(2:7,:));
SDa(3,:) = sum(N(8:15,:));
SDa(4,:) = sum(N(16:21,:));
SDa(5,:) = sum(N(22:end,:));

% scale to create size distribution:
SD = SD./repmat(sum(SD),[size(L,1),1]);
SDa = SDa./repmat(sum(SDa),[size(L,1),1]);


figure(1)
clf
set(gcf,'units','cent','position',[10    7.9022    21   20])
hold on
PlotStage = [1:23 52];
SPs = reshape(1:24,[3,8])';
for i = 1:length(PlotStage)
%subplot(length(PlotStage),1,i)
subplot(8,3,SPs(i))
Tmp = [SD(:,PlotStage(i)), SDa(:,PlotStage(i))];
bar(Tmp,1,'stacked')
text((1:5)-0.2,SD(:,PlotStage(i))+0.1,num2str(round(SD(:,PlotStage(i)),3)),'fontsize',8)
text((1:5)+0.2,SDa(:,PlotStage(i))+0.1,num2str(round(SDa(:,PlotStage(i)),3)),'fontsize',8)
set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'xtick',1:5,'xticklabel',[],'ytick',0:0.25:1,'ylim',[0 1.2])
end
xlabel('Stage','fontsize',14)
ylabel('Proportion of cohort','fontsize',14)


%37% mortality reduction
Surv2 = Surv;
Surv2(9:end) = 1 - (1-Surv(9:end))*(1-0.37);
LA2 = LA;
LA2(2:end,1:end-1) = diag(Surv2);

% Start off at SAD
[Wa,La] = eig(LA);
La = diag(La);
Wa = Wa(:,La==max(La));

N = nan(length(Wa),50);
N(:,1) = Wa./Wa(1).*1e6;

for t = 2:50
    if t < 7
    N(:,t) = LA*N(:,t-1);
    else
    N(:,t) = LA2*N(:,t-1);
    end
end

figure(2)
subplot(2,1,2)
hold on
plot(1:50,sum(N(end-Ts(end):end,:)),'k-')

La = Lmaster;
Survs = [0.6747, 0.75, 0.6758, 0.7425, 0.8091]; % from Table 1 in Crowder
La(3:4,3) = La(3:4,3)./Survs(3).*(1-(1-Survs(3)).*0.37);
La(4:5,4) = La(4:5,4)./Survs(4).*(1-(1-Survs(4)).*0.37);
La(5,5) = La(5,5)./Survs(5).*(1-(1-Survs(5)).*0.37);

Ns = nan(5,50);
Ns(:,1) = W./W(1).*1e6;

for t = 2:50
    if t < 7
    Ns(:,t) = Lmaster*Ns(:,t-1);
    else
    Ns(:,t) = La*Ns(:,t-1);
    end
end

plot(1:50,Ns(5,:),'k--')

ylim([0 1.5e5])

set(gca,'tickdir','out','ticklength',[0.01 0.01])
set(gca,'ylim',[0 1.5e5],'ytick',0:5e4:1.5e5)
set(gca,'fontsize',10)
xlabel('Time (years)','fontsize',12)
ylabel('Number of adult female turtles','fontsize',12)

        
% Also add in actual sizes to size distribution?

% Elasticities of age-model
% first get the dominant eigenvectors:
[W,Ls,V]=eig(LA);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
V = V(:,maxL);
l = max(Ls);

sig = [Surv 0];
sig = sig(:);
% Stage durations:

Ts = 1;
l = max(Ls);
% sensitivity of growth to survival:
sGsig = (sig./l).^(Ts-1)./(( (sig./l).^Ts - 1).^2) .* (Ts.*(1-sig./l) +  (sig./l).^(Ts+1) -(sig./l));
sGsig(end) = 0;
% sensitivity of within-stage-survival to survival
sPsig = 1 - sGsig;
% survival-within-stage sensitivity
sP = V(:).*W(:)./dot(V,W);
% growth senstivity
sG = V(2:end).*W(1:end-1)./dot(V,W);
sG(end+1) = 0;

% elasticity to sigma:
E = sig(:)./l.*(sP(:).*sPsig(:) + sG(:).*sGsig(:));

% Combine elasticities by stage:
%Ts = [1 7 8 6 32];
Combo_mat = zeros(5,length(sig));
Combo_mat(1,1) = 1;
Combo_mat(2,2:8)= 1;
Combo_mat(3,9:16)=1;
Combo_mat(4,17:(17+6))=1;
Combo_mat(5,24:end) = 1;

EE = Combo_mat*E(:);

figure; 
subplot(2,1,1)
bar(E)
subplot(2,1,2)
bar(EE)




