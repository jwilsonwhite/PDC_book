function Dungeness

% Model to look at Dungeness crab dynamics, including nemertean worm effect
% Based on model of Hobbs & Botsford (1989) TPB

doLs = false; % switch to do calculations of eigenvalues, or not (model with worm only)
doWCs = false; % switch to do calculations of eigenvalues (modwl with worm + cannibalism)

A = 1:13; % Age classes
% Influence functions, from Fig. 3 in Hobbs & Botsford (1989)
Infl_b = [0 0 0.24 0.26 0.28 0.12 0.05 0.02 0.01 0.005 0.005 0 0];
% Normalize:
Infl_b = Infl_b./sum(Infl_b);

% Influence on cannibalism, from Botsford & Hobbs (1995)
Infl_c = [0 0.025 0.18 0.325 0.3 0.12 0.05 0.02 0.01 0 0 0 0];
% Normalize
Infl_c = Infl_c./sum(Infl_c);

% values of K and G to vary:
Ks =  1:-0.05:-3;
Gs = 0.025:0.025:0.975;

% Eigenvalue calculation for stability:
if doLs
syms L
Lans = nan(length(Ks),length(Gs));
for k = 1:length(Ks)
    for g = 1:length(Gs)
        
        % The worm influence function is scaled by overall lifetime egg
        % production (Phi_b)
        % The survival at equilibrium (G) is 1/Phi_b.
        Infl_w = [NaN, Infl_b(1:end-1)/Gs(g)]./(1+1/Gs(g));
        Infl_w(1) = -1/(1+1/Gs(g));
        
        % characteristic equation:
        X = sum( (Infl_b + Ks(k)*Infl_w).*L.^-A);
        S  =solve(X==1,L,'ReturnConditions',true);
        LL = vpa(S.L); % numerically evaluate the roots
        LLind = find((abs(LL)==max(abs(LL))),1); % find the dominant eigenvalue
        if isempty(LLind)
            Lans(k,g) = NaN;
        else
        Lans(k,g) = LL(LLind);
        end
        
    end % end g loop
end % end k loop
save Dungeness_worm_Ls.mat Lans Ks Gs
else
load Dungeness_worm_Ls.mat Lans Ks Gs
end

% Plot stability boundaries:
Mag = abs(Lans) < 1;
Theta  = acos(real(Lans)./abs(Lans)) == 0;
Stab = nan(size(Lans));
Stab(~Mag & Theta) = 1;
Stab(Mag & Theta) = 2;
Stab(Mag & ~Theta) = 3;
Stab(~Mag & ~Theta) = 4;

% Figure:
figure(1)
clf
set(gcf,'units','cent','position',[10 10 18 12])

% subplot 1: eigenvalue magnitude
subplot(2,2,1)
hold on
colormap(1-gray(8))
contourf(Gs,Ks,abs(Lans),-1:0.1:3) % magnitude of dominant eigenvalue
%contour(Gs,Ks,Stab,3) % stability regions
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[-3 1],'ytick',-3:1:1,'xlim',[0,1],'xtick',0:0.2:1)
set(gca,'fontsize',10)
xlabel('Equilibrium egg survival rate (g[W*])','fontsize',12)
ylabel('Elasticity of egg predation rate (K_w)','fontsize',12)
ch = colorbar;
set(ch,'tickdir','out','ticklength',[0.015])
ylabel(ch,'|lambda|')

% subplot 2: period of oscillation
subplot(2,2,3)
hold on
colormap(1-gray)
Th = 2*pi./abs(asin(imag(Lans)./abs(Lans)));
Th(Th>20)=20;
Th(imag(Lans)==0) = NaN;
contourf(Gs,Ks,Th) % period of oscillation (2pi/theta)
%contour(Gs,Ks,Stab,3) % stability regions
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[-3 1],'ytick',-3:1:1,'xlim',[0,1],'xtick',0:0.2:1)
set(gca,'fontsize',10)
xlabel('Equilibrium egg survival rate (g[W*])','fontsize',12)
ylabel('Elasticity of egg predation rate (K_w)','fontsize',12)
ch = colorbar;
set(ch,'tickdir','out','ticklength',[0.015])
ylabel(ch,'Period (y)')

% Relationship between egg survival & worm density (from Fig 6 in H&B 1989)
Worm_dens = [0.05, 0.065, 0.09, 0.235];
Egg_surv = [0.78, 0.72, 0.575, 0.56];

subplot(2,2,2)
hold on
plot(Worm_dens,Egg_surv,'ko')
b = regress(Egg_surv(:),[ones(length(Egg_surv),1),Worm_dens(:)]);
plot(Worm_dens,b(1)+b(2)*Worm_dens,'k-')
xlabel('Normalized worm density','fontsize',12)
ylabel('Egg survival','fontsize',12)
set(gca,'xlim',[0 0.25],'ylim',[0.4 0.9])
set(gca,'ytick',0:0.2:1,'xtick',0:0.1:1)
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'fontsize',10)
xlabel('Normalized worm population size','fontsize',12)
ylabel('Survival rate','fontsize',12)


% Analysis including both worm & cannibalism:
% values of K and G to vary:
Kws =  3:-0.05:-3;
Kcs = -4:0.05:1;

if doWCs
syms L

G = 0.7;
 % The worm influence function is scaled by overall lifetime egg
        % production (Phi_b)
        % The survival at equilibrium (G) is 1/Phi_b.
Infl_w = [NaN, Infl_b(1:end-1)/G]./(1+1/G);
Infl_w(1) = -1/(1+1/G);

Infl_c(2) = 0.02;


Lans = nan(length(Kws),length(Kcs));
for k = 1:length(Kws)
    for kk = 1:length(Kcs)
       
        % characteristic equation:
        X = sum( (Infl_b + Kws(k)*Infl_w + Kcs(kk)*Infl_c).*L.^-A);
        S  =solve(X==1,L,'ReturnConditions',true);
        LL = vpa(S.L); % numerically evaluate the roots
        LLind = find((abs(LL)==max(abs(LL))),1); % find the dominant eigenvalue
        if isempty(LLind)
            Lans(k,kk) = NaN;
        else
        Lans(k,kk) = LL(LLind);
        end
    end % end g loop
end % end k loop
save Dungeness_worm+cann_Ls.mat Lans Kws Kcs
else
load Dungeness_worm+cann_Ls.mat Lans Kws Kcs
end % end if doWCs


% Plot stability boundaries:
Mag = abs(Lans) < 1;
Theta  = asin(imag(Lans)./abs(Lans)) == 0;
Stab = nan(size(Lans));
Stab(~Mag & Theta) = 1;
Stab(Mag & Theta) = 2;
Stab(Mag & ~Theta) = 3;
Stab(~Mag & ~Theta) = 4;

% Conditions for proper oscillations
Th = abs(asin(imag(Lans)./abs(Lans)));
Th(Th>0) = 2*pi./Th(Th>0);

OK = abs(Lans) > 0.95 & Th <= 12 & Th >= 9;


% subplot 4: tongue of oscillations
subplot(2,2,4)
hold on
colormap(1-gray)
contourf(Kcs,Kws,OK,1) % period of oscillation (2pi/theta)
contour(Kcs,Kws,Stab,3) % stability regions
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(gca,'ylim',[-3 3],'ytick',-3:1:3,'xlim',[-4,1],'xtick',-4:1:1)
set(gca,'fontsize',10)
xlabel('Elasticity of cannibalism rate (K_c)','fontsize',12)
ylabel('Elasticity of egg predation rate (K_w)','fontsize',12)
plot([-4 1],[0 0],'k--')
plot([-4 1],[-1.1 -1.1],'k--')


