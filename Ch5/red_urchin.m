function red_urchin

% Simulate population model used by Botsford et al. (1994) for red sea
% urchin
% Illustrates how multi-modal size distributions can arise without fancy
% mechanisms

% Parameters for red sea urchin:
k = 0.3; % vB growth
Linf = 100; % vB asymptotic maximum
sig_k = 0.01; % SD of k
sig_Linf = 100; % SD of Linf

A = 0:50; % age classes
L = 0:0.1:99; % length classes
Rmu = 60; % mean recruitment

%------------------------------------
% Plot results (no variability case):
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 18])


% Solution to size-based von Foerster equation:
mu = 0;
sd = NaN;
p = [k, Linf, mu, sd, Rmu];
N = SizeStruct_VarR_EqVB(1,p,A,L,'case1');

%%%%%%%%%%%%%
subplot(2,1,1)
hold on

plot(L,N,'k-','linewidth',1);
ylim([0 200])
xlim([0 100])

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10)
ylim([0 200])
set(gca,'xtick',0:20:100,'xlim',[0 100])
ylabel('Abundance (n[m])','fontsize',12)
xlabel('Test diameter (mm)','fontsize',12)
%%%%%%%%%%%%%

%%%%%%%%%%%%%
subplot(2,1,2)
hold on

plot(L,N,'k-','linewidth',1);

% Pulsed recruitment
mu = 0.5;
sd = 0.15;
p = [k, Linf, mu, sd, Rmu];
N = SizeStruct_VarR_EqVB(1,p,A,L,'case1');

plot(L,N,'k-','color',[0.5 0.5 0.5],'linewidth',1)

ylim([0 200])
xlim([0 100])
set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10)
ylim([0 200])
set(gca,'xtick',0:20:100,'xlim',[0 100])
ylabel('Abundance (n[m])','fontsize',12)
xlabel('Test diameter (mm)','fontsize',12)
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%------------------------------------

% Now, variable recruitment

%------------------------------------
% plot the different growth curves:
Nc = 1e2; %1e5; % (Number of simulations: more = smoother curve)
gamk1 = k^2/sig_k; % shape param
gamk2 = k/gamk1; % scale param 
ks = gamrnd(gamk1,gamk2,[Nc,1]); % distribution of k parameters
gamL1 = Linf^2/sig_Linf; % shape param 
gamL2 = Linf/gamL1; % scale param
Linfs = gamrnd(gamL1,gamL2,[Nc,1]); % distribution of Linf parameters

% Plotting
figure(2)
clf
set(gcf,'units','cent','position',[10 10 9 25])

sh(1) = subplot(3,1,1);
hold on
for i = 1:Nc
    Ls(i,:) = Linfs(i).*(1-exp(-k.*A));
end
plot(A,mean(Ls),'k-')
plot(A,quantile(Ls,0.025),'k-')
plot(A,quantile(Ls,0.975),'k-')

sh(2) = subplot(3,1,2);
hold on
for i = 1:Nc
    Ls(i,:) = Linf.*(1-exp(-ks(i).*A));
end
plot(A,mean(Ls),'k-')
plot(A,quantile(Ls,0.025),'k-')
plot(A,quantile(Ls,0.975),'k-')

sh(3) = subplot(3,1,3);
hold on
for i = 1:Nc
    Ls(i,:) = Linfs(i).*(1-exp(-ks(i).*A));
end
plot(A,mean(Ls),'k-')
plot(A,quantile(Ls,0.025),'k-')
plot(A,quantile(Ls,0.975),'k-')

set(sh(:),'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10,...
    'ylim',[0 150],'xtick',0:10:30,'xlim',[0 30])
ylabel(sh(2),'Test diameter (mm)','fontsize',12)
xlabel(sh(3),'Age (y)','fontsize',12)

%------------------------------------

%------------------------------------
% Figure 3: show effect on variable growth on size distribution:
figure(3)
clf
set(gcf,'units','cent','position',[10 10 9 25])

% Some different means & sds of 
mu1 = 0; mu2 = 0.5; % timing of recruitment pulse (as a fraction of year)
sd1 = NaN; sd2 = 0.15; % sd of recruitment pulse
A  = 0:50; % age classes
L = 0:0.1:140; % length classes
% pre-allocate variables
N1 = nan(Nc,length(L));
N2 = N1; N3 = N1; N4 = N1; N5 = N1; N6 = N1; N7 = N1;

for i = 1:Nc
    % constant recruitment runs:
p = [k, Linfs(i), mu1, sd1, Rmu];
N1(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
p = [ks(i), Linf, mu1, sd1, Rmu];
N2(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
p = [ks(i), Linfs(i), mu1, sd1, Rmu];
N3(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');

% periodic recruitment runs:
p = [k, Linfs(i), mu2, sd2, Rmu];
N4(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
p = [ks(i), mean(Linfs), mu2, sd2, Rmu];
N6(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
p = [ks(i), Linfs(i), mu2, sd2, Rmu];
N7(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2');

end
p = [k, Linf, mu1, sd1, Rmu];
N5 = SizeStruct_VarR_EqVB(1,p,A,L,'case2');
p = [k*1.001, Linf*1.001, mu2, sd2, Rmu];
N5b = SizeStruct_VarR_EqVB(1,p,A,L,'case2');

N2(N2==0) = NaN;
N5(N5==0) = NaN;
N6(N6==0) = NaN;

sh(1) = subplot(3,1,1);
hold on

plot(L,nanmean(N1)./nansum(nanmean(N1)),'k','color',[0.5 0.5 0.5],'linewidth',1) % var Linf
plot(L,nanmean(N2)./nansum(nanmean(N2)),'k--','color',[0.5 0.5 0.5],'linewidth',1) % var k
plot(L,nanmean(N3)./nansum(nanmean(N3)),'k-','linewidth',1) % var both
plot(L,N5./nansum(N5),'k-','linewidth',3) % deterministic


sh(2) = subplot(3,1,2);
hold on

plot(L,N4(1,:)./nansum(N4(1,:)),'k-','linewidth',3)
plot(L,nanmean(N4./nansum(nanmean(N4))),'k','color',[0.5 0.5 0.5],'linewidth',1)
plot(L,nanmean(N6)./nansum(nanmean(N6)),'k--','color',[0.5 0.5 0.5],'linewidth',1)
plot(L,nanmean(N7)./nansum(nanmean(N7)),'k-','color','k','linewidth',1)

%------------------------------------

%------------------------------------
% Panel 3
% Variation in mortality schedule

mu1 = 0; 
sd1 = NaN;
A  = 0:50;
L = 0:0.1:140;
N = nan(Nc,length(L));
for i = 1:Nc
p = [ks(i), Linfs(i), mu1, sd1, Rmu];
N1(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case2'); % constant mortality
N2(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case3'); % Lower refuge
N3(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case4'); % upper refuge
N4(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case5'); % lower & upper refuge

p = [ks(i), Linfs(i), 0.5, 0.15, Rmu*2];
N5(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case4'); % upper refuge + pulsed
N6(i,:) = SizeStruct_VarR_EqVB(1,p,A,L,'case5'); % lower & upper refuge + pulsed

end

% 
figure
hold on

plot(L,nanmean(N1)./sum(nanmean(N1)),'color','b')

plot(L,nanmean(N3)./sum(nanmean(N3)),'k--','color','r')

set(sh(:),'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k','fontsize',10,...
    'ylim',[0, 0.004],'xtick',0:20:150,'xlim',[0, 130])
xlabel(sh(3),'Test diameter (mm)','fontsize',12)
ylabel(sh(2),'Abundance (n[m])','fontsize',12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N = SizeStruct_VarR_EqVB(t,p,A,L,mortcase)
% size-based von Foerster equation
% Pulsed recruitment, no variability in von Bertalanffy equation

% read in parameters: 
k = p(1); % vB growth 
Linf = p(2); % vB asymptotic maximum length
mu = p(3); % mean day of year (as a fraction of 365 days) on which recruitment happens
sd = p(4); % sd of recruitment timing
Rmu = p(5); % mean recruit density

% von Bertalanffy growth
L0 = 0;
A=log(1-L/Linf)./-k; % inverse of vB equation
A(L>Linf) = 0;
A(A==0) = max(A); % substitute the maximum value


gl0 = k*(Linf-L0);
gl =  k*(Linf-L);

if isnan(sd) % is recruitment pulsed or constant?
% Constant recruitment
R = Rmu;
else
% Pulsed recruitment. Mean = Jun 30. Time is fractions of a year
tt = mod(t-A,1); % take fraction of the year for all ages going into the past
R = normpdf(tt,mu,sd);
dA = diff(A);
dA = [dA(1) dA]; % pad with a dummy first entry to extend length
R = R./mean(R);
R = R.*Rmu;
end

% Mortality rate: (several different options)
switch mortcase
    case 'case1'
% unimodal mortality rate:
D = 0.60 - 0.59/50^2*(L-50).^2; % min at 0.01, peak of 0.6 at 50 mm 
D = max(D,0.01);
    case 'case2'
% constant mortality rate:
D = 0.1;
D = 0.2;
    case 'case3'
% low mortality for small individuals:
D = ones(size(L)).*0.1;
D(L<40) = 0;
    case 'case4'
% high mortality for small individuals:
%D = ones(size(L)).*0.1;
%D(L<80) = 1;

D = ones(size(L))*0.01;
D(L<80) = 0.2;

    case 'case5'
% high mortality for medium individuals:
D = ones(size(L)).*1.0;
D(L<40) = 0;
D(L>80) = 0.1;
end

% Integrate D/g(l) wrt l
Int = cumsum(D./gl).*diff(L(1:2));

% Solution:
N = R.*gl0./gl.*exp(-Int);
%N(L>0.99*Linf)=0;
N(N<0) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dn = SizeStruct_ConsR_EqVB(l,n,p)
% size-based von Foerster equation
% Constant recruitment, no variability in von Bertalanffy equation

k = p(1);
Linf = p(2);

% von Bertalanffy growth
g = k*(Linf-l);
gl = -k;% ?g/?l = -k

% unimodal mortality rate:
D = 0.60 - 0.59/50^2*(l-50).^2; % min at 0.01, peak of 0.6 at 50 mm 
D = max(D,0.01);

dn = -n*(gl + D)./g;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
