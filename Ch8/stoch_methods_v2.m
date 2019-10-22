function stoch_methods_v2(doSims)

% Compare different methods for estimating the stochastic growth rate
% using Totoaba as an example to parallel what Fieburg & Ellner (2001) did
% This version - check for agreement with extinction probabilities, not
% estimates of mu

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

% Indices:
Find = A*(0:(A-1))+1;
Jind = [2 A+3];
SAind = A.*(2:5)+(4:7);
Aind = A.*(6:23)+(8:25);


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

% get the dominant eigenvectors:
[W,Ls,V]=eig(L);
Ls = diag(Ls);
maxL = Ls==max(Ls);
W = W(:,maxL);
V = V(:,maxL);

N0 = W./sum(W);
%--------------------------------------


%--------------------------------------
% setup simulations
n = 1e3; % number of simulated actual trajectories
CVs = [0.1,0.5]; % CVs of variation
Tmaster = 1e3; % number of actual simulated years to estimate true p(E)
Ts = [2,5,10]; % number of years sampled to estimate L
X = 0.1; % extinction threshold
Thresh = 100; % extinction time threshold

if ~exist('doSims','var')
    doSims=true;
end

% initial run to get true mu
if doSims


% Apply different estimation methods:
for c = 1:length(CVs)
    
for nn = 1:n
    
    
    % Find true extinction rate by simulation
    N0 = W./sum(W);
    N = zeros(A,Thresh);
    N(:,1) = N0;
    
    for t = 2:Thresh
    % variation in L matrix
    Lt = L; 
    Lt(Find) = max(eps,Lt(Find) + Lt(Find).*normrnd(0,CVs(c)));
    Lt(Jind) = min(1,max(eps,Lt(Jind) + Lt(Jind).*normrnd(0,CVs(c))));
    Lt(SAind) = min(1,max(eps,Lt(SAind) + Lt(SAind).*normrnd(0,CVs(c))));
    Lt(Aind) = min(1,max(eps,Lt(Aind) + Lt(Aind).*normrnd(0,CVs(c))));
    
    % store L
    Lts(:,:,t) = Lt;
    
    % Advance model
    N(:,t) = Lt*N(:,t-1);
    
    end % end loop over t
    
% Estimate true mu
%NN = log(sum(N));
%b = regress(NN(:),[ones(Tmaster,1),(1:Tmaster)']);
%mu_true(c) = b(2);

% Estimate actual extinction rate:
Actual_extinct(c,n) = any(sum(N)<X); 

end % end loop over n
    
    
    % Now do estimation:
    % different numbers of samples
    for t = 1:length(Ts)
        
       % RTM
       % grab Ts(t) matrices
       Rand_RTM = randsample(n,Ts(t));
       RTM_mats = Lts(:,:,Rand_RTM);
       
       % PMM
     %  keyboard
       Fecs = squeeze(RTM_mats(1,:,:))';
       Mean_Fec = mean(Fecs);
       Survs = squeeze(sum(RTM_mats(2:end,1:end-1,:)))';
       Mean_Surv = mean(Survs);
       All_rates = [Fecs, Survs];
       Mean_all = mean(All_rates);
       Covmat = cov(All_rates);
     
       % SFA
       Lmean = mean(RTM_mats,3);
       [W,Ls,V]=eig(Lmean);
       Ls = diag(Ls);
       maxL = Ls==max(Ls);
        W = W(:,maxL);
        V = V(:,maxL);
        l = max(Ls);
        % fecundity elasticities:
        eF = V(1).*W'.*Mean_Fec(:)'./l./dot(V,W);
        % survival elasticities:
        eS = V(2:end).*W(1:end-1).*Mean_Surv(:)./l./dot(V,W);
        E = [eF(:);eS(:)];
       
        % SFA calculation
        D = Covmat./repmat(Mean_all(:),[1,length(E)])./repmat(Mean_all(:)',[length(E),1]);
        D(isnan(D)) = 0;
        Sig2_sfa(c,nn,t) = E(:)'*D*E(:);
        mu_sfa(c,nn,t) = log(l) - Sig2_sfa(c,nn,t)/2;
        
        Px_sfa(cv,t) = normcdf( (-log(1/X) - Mu_sfa(nn,t).*Thresh)./sqrt(Sig2(nn,t).*Thresh))...
             + exp((-2.*Mu_sfa(nn,t).*log(1/X))./Sig2_sfa(nn,t)).*...
             (1-normcdf( (log(1/X) - Mu_sfa(nn,t).*Thresh)./sqrt(Sig2_sfa(nn,t).*Thresh)));
        
        % Iterate for RTM and PMM calculations
        N_rtm = zeros(A,Tmaster);
        N_rtm(:,1) = W./sum(W);
        N_pmm = N_rtm;
        for tt = 2:Thresh
            N_rtm(:,tt) = Lts(:,:,randsample(Rand_RTM,1))*N_rtm(:,tt-1);
            
            PMM_vec = mvnrnd(Mean_all,Covmat);
            L_pmm = diag(min(1,max(0,PMM_vec(A+1:end))));
            L_pmm = [L_pmm,zeros(A-1,1)];
            L_pmm = [max(0,PMM_vec(1:A)); L_pmm];
            
            N_pmm(:,tt) = L_pmm*N_pmm(:,tt-1);
            
        end %end loop over time simulation
       
       % NN = log(sum(N_rtm));
       % b = regress(NN(:),[ones(Tmaster,1),(1:Tmaster)']);
       % if ~isreal(b); keyboard; end
       % mu_rtm(c,nn,t) = b(2);
      
        
       % NN = log(sum(N_pmm));
       % b = regress(NN(:),[ones(Tmaster,1),(1:Tmaster)']);
       % if ~isreal(b); keyboard; end
       % if isinf(b); keyboard; end
       % mu_pmm(c,nn,t) = b(2);
        
    end % end loop over Ts
    
%end % end loop over nn
end % end loop over CVs
save stoch_methods_sims.mat
else
load stoch_methods_sims.mat
end

% Summarize results with kernel density functions.
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 18])

x = linspace(-1,1,100);

Colors = repmat([0 0.4 0.6]',[1 3]);
Lty = {'-','--','-.'};
SPs = reshape(1:6,[2,3])'; % subplot addresses

for c = 1:length(CVs)
for t = 1:length(Ts)
subplot(3,2,SPs(t,c))
hold on


    sfa = fitdist(squeeze(mu_sfa(c,:,t))'-mu_true(c),'kernel');
    sfa_pdf = pdf(sfa,x);
    rtm = fitdist(squeeze(mu_rtm(c,:,t))'-mu_true(c),'kernel');
    rtm_pdf = pdf(rtm,x);
   % keyboard
    pmm = fitdist(squeeze(mu_pmm(c,:,t))'-mu_true(c),'kernel');
    pmm_pdf = pdf(pmm,x);
    

    plot(x,sfa_pdf,'linestyle',Lty{1},'color',Colors(1,:))
    plot(x,rtm_pdf,'linestyle',Lty{2},'color',Colors(1,:))
    plot(x,pmm_pdf,'linestyle',Lty{3},'color',Colors(1,:))
    

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
         'xcolor',zeros(3,1),'ycolor',zeros(3,1),...
         'fontsize',10,'xtick',-1:0.5:1,'ytick',[])
     xlabel(gca,'Rate of decline in mean N (mu)','fontsize',12)
     ylabel(gca,'Probability density','fontsize',12)

end
end


% notes to self:

% in text, show results of this simulation. basically the method doesn't
% matter, but SFA leads to learning how things work, RTM is easy and
% requires fewer assumptions and thus is used more. 

% F&E say to use PMM or SFA because you can check sensitivity to
% covariation, and this is key. 

% Update from Metcalf uses IPMs (this takes care of critique Loo has of F&E
% results?). Says to use RTM, because less worry about having to figure out
% the proper distributions to use, and don't have to actually estimate
% covariances? 

% Link back to Ch 1 - tension between realism (RTM? PMM?) and analytical
% understanding (SFA). But luckily SFA seems to work ok anyway? but does
% require more work.

% One thing to point out is that in contrast to Cisneros-Mata, the SFA
% works really well when you use the proper distribution to simulate
% things! If you are brute forcing it with normal distributions, then you
% introduce errors that have to be corrected.


% UPDATE 5/4: Looking at mu is not very good, should look at actual p(E)
% instead. So need to calculate n estimates of that, at a few time points, 
% compare vs the various methods. Probably does not make a big difference!


