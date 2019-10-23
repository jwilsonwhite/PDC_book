function totoaba_ext

% Simulate extinction probabilities for Totoaba

%----------------------------------------------
% Leslie matrix
A = 25; % 25 year classes
S = [0.01,0.01 repmat(0.764,1,22)]; % length A-1 survival vector  (0.798?)
Linf = 169.9; k = 0.152; t0 = -0.61;
Len = Linf*(1 - exp(-k*((1:A)-t0))); % length-at-age, from Cisneros-Mata et al. (1995b, Con Bio)
F = Len.^3; % assume fecundity proportional to length
F(1:6) = 0; % immature ages

% Assemble Leslie matrix:
L = [diag(S), zeros(length(S),1)];
L = [F; L];

% Rescale fecundity so that lambda = 1
l = max(eig(L)); % lambda
d = abs( l - 1); % difference between lambda & 1

while d > 1e-4
    if l > 1
        L(1,:) = L(1,:)*0.99; % downscale
    else
        L(1,:) = L(1,:)*1.01; % downscale
    end
    l = max(eig(L)); % lambda
    d = abs( l - 1); % difference between lambda & 1
end


%--------------------------------------
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
%--------------------------------------


%----------------------------------------------

%----------------------------------------------
% Simulations:
X = 0.5; % extinction threshold
T = 100; % time horizon
CV = 0.5:0.5:2;
n = 1e1; % number of simulations

% Set up the various random distributions:
% Juvenile survival:
Ind_juv = [2 A+3]; 
Mean_juv = L(Ind_juv(1));

% Subadult survival:
Ind_sa = A.*(2:5)+(4:7); 
Mean_sa = L(Ind_sa(1));

% Adult survival:
Ind_a = A.*(6:(A-2))+(8:A); 
Mean_a = L(Ind_a(1));

% Reproduction:
Ind_r = A.*(0:A-1)+1;
Mean_r = L(Ind_r);

% Random variability in the parameters, using a normal distribution:
for c = 1:length(CV)
Rands = normrnd(0,repmat(CV(c),[1,T*n*5])); 

% juvenile
Rands_juv = Mean_juv + Mean_juv*Rands;
OK = Rands_juv >= 0 & Rands_juv <=1;
Rands_juv = Rands_juv(OK); % apply constraints
RandsC_juv(c,1:T,1:n) = reshape(Rands_juv(1:(T*n)),[1,T,n]); % only take the first n entries

% subadult
Rands_sa = Mean_sa + Mean_sa*Rands;
OK = Rands_sa >=0 & Rands_sa <=1;
Rands_sa = Rands_sa(OK); % apply constraints
RandsC_sa(c,1:T,1:n) = reshape(Rands_sa(1:(T*n)),[1,T,n]); % only take the first n entries

% Adult
Rands_a = Mean_a + Mean_a*Rands;
OK = Rands_a >= 0 & Rands_a <= 1;
Rands_a = Rands_a(OK); % apply constraints
RandsC_a(c,1:T,1:n) = reshape(Rands_a(1:(T*n)),[1,T,n]); % only take the first n entries


% Repro:
for a = 1:A
Rands_r = Mean_r(a)+Mean_r(a)*Rands;
OK = Rands_r >= 0 ;
Rands_r = Rands_r(OK); % apply constraints
RandsC_r(c,1:T,1:n,a) = reshape(Rands_r(1:(T*n)),[1,T,n]); % only take the first n entries
end
end

% Random variability, using a proper distribution (beta for survival)
% juvenile:
Var = (CV.*Mean_juv).^2;
Alfa = -(Mean_juv^3 - Mean_juv^2 + Var*Mean_juv)./Var;
Alfa = max(realmin,Alfa);
Beta = Alfa*(1/Mean_juv-1);
Rands2_juv = betarnd(repmat(Alfa(:),[1,T,n]),repmat(Beta(:),[1,T,n]));

% subadult:
Var = (CV.*Mean_sa).^2;
Alfa = -(Mean_sa^3 - Mean_sa^2 + Var*Mean_sa)./Var;
Alfa = max(realmin,Alfa);
Beta = Alfa*(1/Mean_sa-1);
Rands2_sa = betarnd(repmat(Alfa(:),[1,T,n]),repmat(Beta(:),[1,T,n]));

% Adult:
Var = (CV.*Mean_a).^2;
Alfa = -(Mean_a^3 - Mean_a^2 + Var*Mean_a)./Var;
Alfa = max(realmin,Alfa);
Beta = Alfa*(1/Mean_a-1);
Rands2_a = betarnd(repmat(Alfa(:),[1,T,n]),repmat(Beta(:),[1,T,n]));


% Reproduction (lognormal)
MeanMeanR = mean(Mean_r);
% not a good analytical expression relating the CV of the lognormal to its
% parameters. So use brute force.
SDs = linspace(0,2,1e3);
Mu = log(MeanMeanR);
CVlist = ( (exp(SDs.^2) - 1).*exp(2.*Mu+SDs.^2) ).^0.5./exp(Mu+SDs.^2/2); 
for c = 1:length(CV)
    CV0 = abs(CVlist-CV(c));
    SD(c) = CVlist(CV0==min(CV0));
end
Rands2_r = normrnd(0,repmat(SD(:),[1,T,n]));

  
% Do simulations
N_juv = zeros(A,T,n);
N_juv(:,1,:) = repmat(W(:),[1,n]); % initial conditions
N2_juv = N_juv; % for use with proper random distribution
N_sa = N_juv; %
N2_sa = N_juv; % for use with proper random distribution
N_a = N_juv; %
N2_a = N_juv; % for use with proper random distribution
N_r = N_juv; %
N2_r = N_juv; % for use with proper random distribution

for c = 1:length(CV)
for nn = 1:n
for t = 2:T
    
    % juvenile:
    Lt = L;
    Lt(Ind_juv) = RandsC_juv(c,t,nn);
    N_juv(:,t,nn) = Lt*N_juv(:,t-1,nn);
    
    Lt = L;
    Lt(Ind_juv) = Rands2_juv(c,t,nn);
    N2_juv(:,t,nn) = Lt*N2_juv(:,t-1,nn);
    
    % subadult:
    Lt = L;
    Lt(Ind_sa) = RandsC_sa(c,t,nn);
    N_sa(:,t,nn) = Lt*N_sa(:,t-1,nn);
    
    Lt = L;
    Lt(Ind_sa) = Rands2_sa(c,t,nn);
    N2_sa(:,t,nn) = Lt*N2_sa(:,t-1,nn);
    
    % adult:
    Lt = L;
    Lt(Ind_a) = RandsC_a(c,t,nn);
    N_a(:,t,nn) = Lt*N_a(:,t-1,nn);
    
    Lt = L;
    Lt(Ind_a) = Rands2_a(c,t,nn);
    N2_a(:,t,nn) = Lt*N2_a(:,t-1,nn);
    
    % repro:
    Lt = L;
    Lt(Ind_r) = RandsC_r(c,t,nn,:);
    N_r(:,t,nn) = Lt*N_r(:,t-1,nn);
    
    Lt = L;
    Ltmp = log(Lt(Ind_r));
    Ltmp = Ltmp + Rands2_r(c,t,nn);
    Ltmp = exp(Ltmp);
    Lt(Ind_r) = Ltmp;
    N2_r(:,t,nn) = Lt*N2_r(:,t-1,nn);
  
end % end loop over T
end % end loop over n

% Quasi-extinctions:
NN = squeeze(sum(N_juv,1));
Ex = sum(NN<=X)>0;
Px_juv(c) = mean(Ex);

NN = squeeze(sum(N2_juv,1));
Ex = sum(NN<=X)>0;
Px2_juv(c) = mean(Ex);

NN = squeeze(sum(N_sa,1));
Ex = sum(NN<=X)>0;
Px_sa(c) = mean(Ex);

NN = squeeze(sum(N2_sa,1));
Ex = sum(NN<=X)>0;
Px2_sa(c) = mean(Ex);

NN = squeeze(sum(N_a,1));
Ex = sum(NN<=X)>0;
Px_a(c) = mean(Ex);

NN = squeeze(sum(N2_a,1));
Ex = sum(NN<=X)>0;
Px2_a(c) = mean(Ex);

NN = squeeze(sum(N_r,1));
Ex = sum(NN<=X)>0;
Px_r(c) = mean(Ex);

NN = squeeze(sum(N2_r,1));
Ex = sum(NN<=X)>0;
Px2_r(c) = mean(Ex);

end % end loop over CV

save totoaba_ext_results.mat -v7.3

% Calculation for analytical solution:
% Juvenile:
Sig2 = (E(26:27)'*ones(2)*E(26:27)).*CV.^2;
Mu = -Sig2/2;

% Calculation for *corrected* analytical solution (use actual
% post-truncation CV)
RandsC2 = reshape(RandsC_juv,[length(CV),T*n]);
CV_act_juv = std(RandsC2,[],2)./mean(RandsC2,2);
H1(1,:) = RandsC2(CV==0.5,:);
% account for effect on mean:
for c = 1:length(CV)
    Lt = L;
    Lt(Ind_juv) = mean(RandsC2(c,:));
    lt(c) = max(eig(Lt));
end
Sig2_act = (E(26:27)'*ones(2)*E(26:27)).*CV_act_juv.^2;
Sig2_act = Sig2_act';
Mu_act = log(lt)-Sig2_act/2;

Px_an_juv = normcdf( (-log(1/X) - Mu.*T)./sqrt(Sig2.*T))...
             + exp((-2.*Mu.*log(1/X))./Sig2).*...
             (1-normcdf( (log(1/X) - Mu.*T)./sqrt(Sig2.*T)));
         
Px_an_act_juv = normcdf( (-log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))...
             + exp((-2.*Mu_act.*log(1/X))./Sig2_act).*...
             (1-normcdf( (log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T)));
         
 % Subadult:
Sig2 = (E(28:31)'*ones(4)*E(28:31)).*CV.^2;
Mu = -Sig2/2;

% Calculation for *corrected* analytical solution (use actual
% post-truncation CV)
RandsC2 = reshape(RandsC_sa,[length(CV),T*n]);
CV_act_sa = std(RandsC2,[],2)./mean(RandsC2,2);
H1(2,:) = RandsC2(CV==0.5,:);
% account for effect on mean:
for c = 1:length(CV)
    Lt = L;
    Lt(Ind_sa) = mean(RandsC2(c,:));
    lt(c) = max(eig(Lt));
end

Sig2_act = (E(28:31)'*ones(4)*E(28:31)).*CV_act_sa.^2;
Sig2_act = Sig2_act';
Mu_act = log(lt)-Sig2_act/2;

Px_an_sa = normcdf( (-log(1/X) - Mu.*T)./sqrt(Sig2.*T))...
             + exp((-2.*Mu.*log(1/X))./Sig2).*...
             (1-normcdf( (log(1/X) - Mu.*T)./sqrt(Sig2.*T)));
         
Px_an_act_sa = normcdf( (-log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))...
             + exp((-2.*Mu_act.*log(1/X))./Sig2_act).*...
             (1-normcdf( (log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T)));        
         
% Adult:
Sig2 = (E(32:end)'*ones(18)*E(32:end)).*CV.^2;
Mu = -Sig2/2;

% Calculation for *corrected* analytical solution (use actual
% post-truncation CV)
RandsC2 = reshape(RandsC_a,[length(CV),T*n]);
CV_act_a = std(RandsC2,[],2)./mean(RandsC2,2);
H1(3,:) = RandsC2(CV==0.5,:);

% account for effect on mean:
for c = 1:length(CV)
    Lt = L;
    Lt(Ind_a) = median(RandsC2(c,:));
    lt(c) = max(eig(Lt));
end

Sig2_act = (E(32:end)'*ones(18)*E(32:end)).*CV_act_a.^2;
Sig2_act = Sig2_act';
Mu_act = log(lt)-Sig2_act/2;

Px_an_a = normcdf( (-log(1/X) - Mu.*T)./sqrt(Sig2.*T))...
             + exp((-2.*Mu.*log(1/X))./Sig2).*...
             (1-normcdf( (log(1/X) - Mu.*T)./sqrt(Sig2.*T)));
         
Px_an_act_a = normcdf( (-log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))...
             + exp((-2.*Mu_act.*log(1/X))./Sig2_act).*...
             (1-normcdf( (log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))); 
          
         
% Repro:
Sig2 = (E(1:25)'*ones(25)*E(1:25)).*CV.^2;
Mu = -Sig2/2;

% Calculation for *corrected* analytical solution (use actual
% post-truncation CV)
RandsC2 = reshape(RandsC_r,[length(CV),T*n,A]);
CV_act_r = std(RandsC2,[],2)./mean(RandsC2,2);

Ht = median(RandsC2,3);
H1(4,:) = Ht(CV==0.5,:);
% account for effect on mean:
for c = 1:length(CV)
    Lt = L;
    Lt(Ind_r) = mean(RandsC2(c,:,:),2);
    lt(c) = max(eig(Lt));
end


Sig2_act = (E(1:25)'*ones(25)*E(1:25)).*squeeze(nanmean(CV_act_r,3)).^2;
Sig2_act = Sig2_act';
Mu_act = log(lt)-Sig2_act/2;

Px_an_r = normcdf( (-log(1/X) - Mu.*T)./sqrt(Sig2.*T))...
             + exp((-2.*Mu.*log(1/X))./Sig2).*...
             (1-normcdf( (log(1/X) - Mu.*T)./sqrt(Sig2.*T)));
         
Px_an_act_r = normcdf( (-log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))...
             + exp((-2.*Mu_act.*log(1/X))./Sig2_act).*...
             (1-normcdf( (log(1/X) - Mu_act.*T)./sqrt(Sig2_act.*T))); 
         
                  
%------------------------------         
figure(1)
clf
set(gcf,'units','cent','position',[10 10 21 27])

sh(4)=subplot(4,3,10);
hold on
plot(CV,Px_r,'k')
plot(CV,Px2_r,'color',[0.5 0.5 0.5])
plot(CV,Px_an_r,'k--');
plot(CV,Px_an_act_r,'k-.');


sh(1)=subplot(4,3,1); 
hold on
plot(CV,Px_sa,'k')
plot(CV,Px2_sa,'color',[0.5 0.5 0.5])
plot(CV,Px_an_sa,'k--');
plot(CV,Px_an_act_sa,'k-.');

sh(2) = subplot(4,3,4);
hold on
plot(CV,Px_juv,'k')
plot(CV,Px2_juv,'color',[0.5 0.5 0.5])
plot(CV,Px_an_juv,'k--');
plot(CV,Px_an_act_juv,'k-.');

sh(3) = subplot(4,3,7);
hold on
plot(CV,Px_a,'k')
plot(CV,Px2_a,'color',[0.5 0.5 0.5])
plot(CV,Px_an_a,'k--');
plot(CV,Px_an_act_a,'k-.');


set(sh(1:4),'tickdir','out','ticklength',[0.015 0.015],...
         'xcolor',zeros(3,1),'ycolor',zeros(3,1),...
         'ylim',[0 1],'ytick',0:0.5:2,'xlim',[0 2],'xtick',0:0.5:2,'fontsize',10)
 ylabel(sh(2),'Probability of quasi-extinction by time t, pt(X|N0)','fontsize',12)
 xlabel(sh(4),'Coefficient of variation (n)','fontsize',12)

 sh(5)=subplot(4,3,2);
 histogram(H1(1,:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(6)=subplot(4,3,5);
 histogram(H1(2,:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(7)=subplot(4,3,8);
 histogram(H1(3,:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(8)=subplot(4,3,11);
 histogram(H1(4,:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 
 sh(9)=subplot(4,3,3);
 H = Rands2_juv(CV==0.5,:,:);
 histogram(H(:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(10)= subplot(4,3,6);
 H = Rands2_sa(CV==0.5,:,:);
 histogram(H(:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(11)= subplot(4,3,9);
 H = Rands2_a(CV==0.5,:,:);
 histogram(H(:),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
 sh(12)= subplot(4,3,12);
 H = Rands2_r(CV==0.5,:,:);
 H = log(mean(L(1,:))) + H;
 histogram(exp(H(:)),100,'facecolor',[0.5 0.5 0.5],'normalization','pdf')
 
set(sh(5:12),'tickdir','out','ticklength',[0.015 0.015],...
         'xcolor',zeros(3,1),'ycolor',zeros(3,1),...
         'fontsize',10,'xtick',[])
set(sh([5,9]),'xlim',[0 0.05],'xtick',0:0.01:0.05)
set(sh([6,7,10,11]),'xlim',[0 1],'xtick',0:0.2:1)
set(sh([8,12]),'xlim',[0 6e4],'xtick',0:1e4:6e4)


%-----------------------------------------------
% Fig 2: simpler
figure(1)
clf
set(gcf,'units','cent','position',[10 10 10 27])

sh(4)=subplot(4,1,4);
hold on
plot(CV,Px2_r,'color',[0.5 0.5 0.5])
plot(CV,Px_an_r,'k--');


sh(1)=subplot(4,1,1); 
hold on
plot(CV,Px2_sa,'color',[0.5 0.5 0.5])
plot(CV,Px_an_sa,'k--');


sh(2) = subplot(4,1,2);
hold on
plot(CV,Px2_juv,'color',[0.5 0.5 0.5])
plot(CV,Px_an_juv,'k--');

sh(3) = subplot(4,1,3);
hold on
plot(CV,Px2_a,'color',[0.5 0.5 0.5])
plot(CV,Px_an_a,'k--');


set(sh(1:4),'tickdir','out','ticklength',[0.015 0.015],...
         'xcolor',zeros(3,1),'ycolor',zeros(3,1),...
         'ylim',[0 1],'ytick',0:0.5:2,'xlim',[0 2],'xtick',0:0.5:2,'fontsize',10)
 ylabel(sh(2),'Probability of quasi-extinction by time t, pt(X|N0)','fontsize',12)
 xlabel(sh(4),'Coefficient of variation (n)','fontsize',12)

 

