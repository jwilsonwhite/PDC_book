function dd_eigenvalues

% Generic model to examine unstable 2T cycles

% Life history parameters for Dungeness crab (Hobbs & Botsford 1989)
A = 1:13; % Age classes
Wid = [40 120 150 180 210 240 repmat(240,[1,7])]; % carapace width of each age class
Wgt = Wid.^2.76; % mass as a function of width
SurvF = [repmat(0.2,1,3) repmat(0.9,1,9)]; % female age-specific survival
SurvM = [repmat(0.2,1,6) repmat(0.9,1,6)]; % male age-specific survival
SurvF = [1, cumprod(exp(-SurvF))]; % cumulative survival (female)
SurvM = [1, cumprod(exp(-SurvM))]; % cumulative survival (male)
Surv = mean([SurvF;SurvM]); % average survival across sexes
Eggs = Wgt; % fecundity is proportional to biomass
Eggs(A<3)=0; % exclude immature age classes

% Influence functions
Infl_c = Surv.*Wgt.^0.8; % Influence function of adult biomass on recruit cannibalism
Infl_b = SurvF.*Eggs; % % Influence function on reproduction
% Normalize:
Infl_b = Infl_b./sum(Infl_b);
Infl_c = Infl_c./sum(Infl_c);

%---------------------------------------------------------------
% FIGURE 1
% Illustrate the influence functions + K
figure(1)
clf
set(gcf,'units','cent','position',[10 10 18 10.5])

% Subplot 1: birth influence function
subplot(2,2,1)%
hold on
% plot components of Infl_b (not shown in book)
%plot(A,SurvF./sum(SurvF),'k-','color',[0.5 0.5 0.5])
%plot(A,Eggs./sum(Eggs),'k--','color',[0.5 0.5 0.5])
plot(A,Infl_b,'k')
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[0 0.4],'ytick',0:0.1:0.4,'xlim',[1,10],'xtick',[1 4 7 10 ])
set(gca,'fontsize',10,'xcolor','k','ycolor','k')
ylabel({'Influence on';'reproduction'},'fontsize',12)

% Subplot2: compensation influence function
subplot(2,2,3)
hold on
% plot components of Infl_c (not shown in book)
%plot(A,Surv./sum(Surv),'k-','color',[0.5 0.5 0.5])
%WgtE = Wgt.^0.8;
%plot(A,WgtE./sum(WgtE),'k--','color',[0.5 0.5 0.5])
plot(A,Infl_c,'k-')
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[0 0.4],'ytick',0:0.1:0.4,'xlim',[1,10],'xtick',[1 4 7 10])
set(gca,'fontsize',10,'xcolor','k','ycolor','k')
ylabel({'Influence on';'compensatory mortality'},'fontsize',12)

% Subplots 3-4; Do for both Ricker & Bev-Holt stock-recruit functions

% Choose plotting values at particular values of LEP
LEP = 1e2;
FLEPs = linspace(0.11,1,4);

% Subplot 3: B-H
a = 0.1; % initial slope
b = 1; % max density
N = 0:0.01:10; % range of starting densities
R = a*N./(1 + b.*N); % B-H function
X = linspace(0,1,length(R)); % abscissa values

% Find the equilibria to plot as vertical lines;
Rs = a/b - 1./(b*FLEPs*LEP);
Ns = Rs.*FLEPs;
K = -b .* Ns.*LEP./(1 + b.*Ns.*LEP);
Ks = -b .* N./(1 + b.*N);
for k = 1:length(K)
    K_ind(k) = find(min(abs(Ks-K(k)))==abs(Ks-K(k)));
end

subplot(2,2,2)
hold on
plot(X,R,'k-') % y-axis = recruitment
plot(X,R./N,'k-','color',[0.5 0.5 0.5]) % y-axis = survival
for k = 1:length(K) % plot vertical lines corresponding to the equilibria at representative values of K
    plot([X(K_ind(k)), X(K_ind(k))],[0, 0.1],'k-','color',[0.5 0.5 0.5],'linewidth',0.5)
text(X(K_ind(k)),0.11,num2str(round(K(k),2)));
end
ylabel({'Recruit density';'or survival'},'fontsize',12)
set(gca,'xtick',[],'ytick',[],'xcolor','k','ycolor','k')

% Subplot 4: Ricker
a = 0.2; % initial slope
b = 1; % DD parameter
N = 0:0.01:10; % range of N values
R = a*N.*exp(-b.*N); % Ricker function

% Find the equilibria to plot as vertical lines;
Rs = log(a*FLEPs.*LEP)./(b.*FLEPs.*LEP);
Ns = Rs.*FLEPs;
K = -b*Ns.*LEP;
Ks = -b.*N;
for k = 1:length(K)
    K_ind(k) = find(min(abs(Ks-K(k)))==abs(Ks-K(k)));
end

% Rescale the functions to plot on similar axes
Surv = R./N;
Surv = Surv/max(Surv)*0.1; % Max value will appear at an appropriate y-axis value

subplot(2,2,4)
hold on
plot(X,R,'k-')
plot(X,Surv,'k-','color',[0.5 0.5 0.5])
for k = 1:length(K) % plot vertical lines corresponding to the equilibria at representative values of K
    plot([X(K_ind(k)), X(K_ind(k))],[0, 0.1],'k-','color',[0.5 0.5 0.5],'linewidth',0.5)
text(X(K_ind(k))-0.05,0.11,num2str(round(K(k),2)));
end
ylabel('Offspring survival','fontsize',12)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xtick',[],'ytick',[],'xlim',[0 0.5])
set(gca,'xcolor','k','ycolor','k')

%---------------------------------------------------------------

%---------------------------------------------------------------
% Find the critical values of K for instability

% In theory we could use the symbolic toolbox, but it is faster to find a
% brute-force solution
Ks = -8:0.5:0; % possible values of K
Ms = 0:0.1:1.2; % Possible values of fishing mortality

% Probably a way to do this in a vectorized way, but it is a lot of
% dimensions...
doRuns = false; % a simple switch in case the runs are complete but we are fiddling with plot options
if doRuns
for k = 1:length(Ks)
    for m = 1:length(Ms)
      
SurvM = [repmat(0.2,1,3),repmat(Ms(m)+0.2,1,3), repmat(Ms(m)+0.9,1,6)];
SurvF = [repmat(0.2,1,3),repmat(0.9,1,3), repmat(0.9,1,6)];

SurvF = [1, cumprod(exp(-SurvF))];
SurvM = [1, cumprod(exp(-SurvM))];
      
   Infl_b = SurvF.*Eggs;
   Infl_c = mean([SurvF;SurvM]).*Wgt.^0.8;
   Infl_b = Infl_b./sum(Infl_b);
   Infl_c = Infl_c./sum(Infl_c);
   
   X = sum( (Infl_b + Ks(k).*Infl_c).*L.^-A);
   S  =solve(X==1,L,'ReturnConditions',true);
   LL = vpa(S.L); % numerically evaluate the roots
   LLind = find((abs(LL)==max(abs(LL))),1); % find the dominant eigenvalue
   Lans(k,m) = LL(LLind);
   
   % Some debugging code...
   if abs(double(Lans(k,m)))>= 0.9 && abs(double(Lans(k,m)))<1
     %  keyboard
   end
   
   % Mean of influence function (remove semicolons if you don't want to suppress output):
   Infl = Infl_b + Ks(k).*Infl_c;
   Infl = abs(Infl)/sum(abs(Infl));
   Gen0(k,m) = Infl*A';

   % Median:
   Gen0Med(k,m)= A(find(cumsum(Infl)>=0.5,1));

   end
end
save dd_eig_runs.mat
else
load dd_eig_runs.mat Lans Gen0 Gen0Med
end
%---------------------------------------------------------------


%---------------------------------------------------------------
% Figure with eigenvalues on complex plane, and contour plot of
% eigenvalues (Figure 7.3 in book)

% Plotting:
figure(2)
clf
set(gcf,'units','cent','position',[10 10 10 20])

% Target values of fishing to grab corresponding results:
target1 = 0.0;
target2 = 0.1; % (not used)
target3 = 1.2;

% Find eigenvalues & M that correspond to desired target Ms
LL = double(Lans);
LL = LL(:,Ms>=target1);
Gen0 = Gen0(:,Ms>=target1);

Ms = Ms(Ms>=target1);

% Narrow range of K values that can be searched
L1 = LL(Ks>-3.9&Ks<-3.7,:);
L2 = LL(Ks>-4.1&Ks<-3.9,:);
L3 = LL(Ks>-4.3&Ks<-4.1,:);

% Find the eigenvalues that correspond to the target levels of fishing
% mortality
L1 = LL(:,Ms==target1);
L2 = LL(:,Ms==target3);
Index = 1:2:(length(L1));
L1 = L1(Index);
L2 = L2(Index);

% Plots
subplot(2,1,2)
hold on

% Plot the eigenvalues
plot(real(L1),abs(imag(L1)),'ko','linestyle','-');
plot(real(L2),abs(imag(L2)),'kv','linestyle','-'); 

% Plot unit circle
Theta = linspace(0,2*pi,100);
[x,y] = pol2cart(Theta,1);
plot(x,y,'k-')
plot([-2 2],[0 0],'k--')
plot([0 0],[-2 2],'k--')
set(gca,'xtick',-1:0.5:1,'ytick',-1:0.5:1)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xcolor','k','ycolor','k')
set(gca,'xlim',[-1.3 1.3],'ylim',[-1.3 1.3])
axis square

% Plot contourf of the eigenvalues
subplot(2,1,1)
hold on
Grad = [0 0.4:0.1:1.2];
contourf(Ms,Ks,abs(LL),Grad)
colormap(1-gray(6))
set(gca,'xtick',0:0.2:2)
set(gca,'ytick',-8:2:0)
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015)
axis square
%---------------------------------------------------------------


%---------------------------------------------------------------
% This code isn't used anymore, it plots the real & imaginary parts of the
% eigenvalues
doImag = false;
if doImag
subplot(3,1,2)
hold on
%pcolor(Ms,Ks,double(abs(imag(Lans)))); shading interp
Grad = linspace(0,1.5,10);
contourf(Ms,Ks,abs(imag(LL)),Grad) %[-3:1:1],
%colormap(1-gray(length(Grad)))
colormap(1-gray(6))
set(gca,'xtick',0.1:0.1:0.5,'xlim',[0.1 0.4])
set(gca,'ytick',-8:2:0,'ylim',[-8 0])
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'ytick',0:0.2:1)
axis square


subplot(3,1,2)
hold on

%pcolor(Ms,Ks,Gen0); shading interp
Grad = 2*(3:0.5:5.5);
Grad = 7:0.5:9.5;
%LL = double(Lans);
Gen0(abs(LL)<1)=NaN;
contourf(Ms,Ks,2*Gen0,Grad)
set(gca,'clim',[7 9.5])
%pcolor(Ms,Ks,2*Gen0); shading flat
%colormap(1-gray(length(Grad)))
%colormap(1-gray(6))
set(gca,'xtick',0.1:0.1:0.5,'xlim',[0.1 0.4])
set(gca,'ytick',-8:2:0,'ylim',[-8 0])
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'ytick',1:10)
axis square

subplot(3,1,3)
hold on
%pcolor(Ms,Ks,Gen0); shading interp
%Grad = 2*(3:0.5:5.5);
%LL = double(Lans);
Period = 2*pi./( acos( real(LL)./abs(LL) ));
Period(abs(LL)<1)=NaN;
%Period = double(abs(imag(Lans)))*(2*pi);
%keyboard
contourf(Ms,Ks,Period,Grad)
%pcolor(Ms,Ks,Period); shading flat
%colormap(1-gray(length(Grad)))
%colormap(1-gray(6))
set(gca,'xtick',0.1:0.1:0.5,'xlim',[0.1 0.4])
set(gca,'ytick',-8:2:0,'ylim',[-8 0])
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'ytick',1:10)
axis square
end % end if do Imag
%---------------------------------------------------------------


%---------------------------------------------------------------
% Plots with example timeseries, under increased fishing (Fig. 7.2)

% Simulate example timeseries
T = 1000; % time for simulations
N = zeros(length(A),T); % pre-allocate variables
Nm = N;
Nf = N;

% Function to find equilibrium
eqF = @(F, x, Fp, c) feval(F, x, Fp) - 1/c;

% Recruit survival function
SRfun = @(S,Fp) max(realmin,(1-normcdf(S,5e6,Fp(2)))*Fp(1)); 

% SRfun parameter values that illustrate a range of steepness (and thus
% cyclic behavior)
As = 6.15e-7; % survival at low density
Bs = [2.5e6 1e6 5e5]; % SD of normcdf function
F = 1.2; % instantaneous harvest rate of males

% subplot addresses for plotting
Sp3 = [2 5 8];
Sp2 = [3 6 9];
Sp1 = [1 4 7];

% Set up plot formats
figure(3)
clf
set(gcf,'units','cent','position',[10 10 15 10]) 

figure(4)
clf
set(gcf,'units','cent','position',[5 10 4.5 10]) 

% Loop over SRfun steepness values to plot examples
for i = 1:length(Bs)

% SRfun parameters, for convenience
a = As(1);
b = Bs(i); 

% This figure will hold the nine-panel plot that is Fig. 7.2
figure(3)
subplot(3,3,Sp1(i))
cla
hold on

% recalculate the influence functions, using the fishing rate
SurvM = [repmat(0.2,1,3),repmat(F+0.2,1,3), repmat(F+0.9,1,6)];
SurvF = [repmat(0.2,1,3), repmat(0.9,1,9)];
SurvF = [1, cumprod(exp(-SurvF))];
SurvM = [1, cumprod(exp(-SurvM))];  
   Infl_b = SurvF.*Eggs;
   Infl_c = mean([SurvF;SurvM]).*Wgt.^0.8;
   LEP = sum(Infl_b);
   
Infl_c_uns = sum(Infl_c); % Need an unscaled, summed value for later simulation

% Scale to sum to one
   Infl_b = Infl_b./sum(Infl_b);
   Infl_c = Infl_c./sum(Infl_c);


% First column of plots will show the SR function:
Rs = linspace(0,1e2,1e3); % dummy values of offspring produced
f = SRfun(Rs*Infl_c_uns,[a,b,1]); % apply SRfunction, given Infl_c

% Plot the recruit survival
plot(Rs.*Infl_c_uns,f,'k-');
% Plot the equilibrium
plot([0 max(Rs.*Infl_c_uns)],[1/LEP 1/LEP],'k--')

% Format plot
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'xcolor','k','ycolor','k')
set(gca,'xtick',[],'ytick',[])
set(gca,'xlim',[0 max(Rs.*Infl_c_uns)],'ylim',[0 a*1.1])

% Calculate K (elasticity of slope at equilibrium)
NEq = fzero(@(x) eqF(SRfun,x,[a,b,1],sum(LEP)),1); % Effective N at Eq (scaled by max LEP)
REq = SRfun(NEq,[a,b,1]); % Recruit survival at Eq
% Calculate local slope of SR curve at equilibrium:
Del = 10; % interval for taking derivative
REq2 = feval(SRfun,NEq+Del,[a,b,1]);
slopeAtEq = (REq2-REq)/(Del);
Rstar = NEq./sum(LEP); % equilibrium
K = slopeAtEq*NEq/REq; % elasticity
disp(strcat('K = ',num2str(K)))

   
% Calculate joint influence function, including K
   Infl_c2 = Wgt.^0.8; % this version will be used in iterating the model
   
   Infl_a = Infl_b + K*Infl_c;
   Infl=Infl_a;
   Infl = abs(Infl)/sum(abs(Infl)); % rescale to sum to 1
   Gen = Infl*A'; % generation time
   disp(strcat('Gen = ',num2str(Gen)));
   
   % calculate median of the function:
   min(A(cumsum(Infl)>=0.5))
   
   % Plot just the joint influence function (not used in book)
   figure(4)
   subplot(3,1,i)
   hold on
   plot(A,abs(Infl_a),'k-')
   set(gca,'ylim',[0 max(abs(Infl_a))*1.2])
   set(gca,'xtick',0:5:15,'xlim',[0 max(A)])
   set(gca,'xcolor','k','ycolor','k')
   set(gca,'tickdir','out','ticklength',[0.015 0.015])
   
   % Plot the two influence functions and the joint (Fig. 7.2, middle
   % panels)
   figure(3)
   
   subplot(3,3,Sp3(i))
   hold on
   plot(A,abs(Infl_b),'k-')
   plot(A,K*(Infl_c),'k-')
   plot(A,(Infl_a),'b-')
   plot([1 max(A)],[0 0],'b-')
   plot([1 1],[0 Infl_a(1)],'b-')
   set(gca,'xtick',0:5:15,'xlim',[0 max(A)])
   set(gca,'xcolor','k','ycolor','k')
   set(gca,'tickdir','out','ticklength',[0.015 0.015])
   if i == 1
   set(gca,'ylim',[-1 0.5],'ytick',-1:1:0.5)
   elseif i == 2
   set(gca,'ylim',[-1.1 0.5],'ytick',-3:3)
   else
   set(gca,'ylim',[-3 1],'ytick',-3:1:3)
   end
   
   
   % Third column of Fig. 7.2 will show model simulated timeseries

   % Recalculate unscaled influence funtions for model iteration:
SurvM2 = [repmat(0.2,1,3),repmat(F+0.2,1,3), repmat(F+0.9,1,6)];
SurvF2 = [repmat(0.2,1,3),repmat(0.9,1,3), repmat(0.9,1,6)];

% Initialize variables at equilibrium, pre-allocate settler (S) and recruit
% (R) variables
Nm(:,1) = SurvM*Rstar; Nf(:,1) = SurvF*Rstar;
N=Nm+Nf;
S = zeros(length(N),1);
R = S;

% Loop for simulations
for t = 2:T
    
% Reproduction    
S(t) = Nf(:,t-1)'*Eggs(:);
% Survival
Nf(2:end,t) = Nf(1:end-1,t-1).*exp(-SurvF2)';
Nm(2:end,t) = Nm(1:end-1,t-1).*exp(-SurvM2)';
meanN = (Nf(:,t-1)+Nm(:,t-1))/2;

% Density-dependent recruit survival
R(t) = SRfun(Infl_c2(:)'*meanN,[a,b])*S(t);

% (the following R(t) values should really be divided by 2, and a should be
% multiplied by 2 to obtain same dynamics)
Nf(1,t) = R(t);
Nm(1,t) =R(t);
N=Nm+Nf;

end % end simulations

%disp(R/S) % # recruit survival
%disp(REq) % recruit survival

% Plot sims in final column of Fig 7.2
figure(3)
subplot(3,3,Sp2(i))
plot(450:500,sum(N(2:end,450:500)),'k-')
xlim([450 500])
set(gca,'xtick',450:10:600,'xticklabels',0:10:600)
set(gca,'xcolor','k','ycolor','k')
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[0 200])
set(gca,'ytick',0:50:200)

end % end loop over B values   