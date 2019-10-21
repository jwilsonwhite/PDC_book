function Clark

% Find solutions for the stability analysis of Clark's (1976) delayed
% recruitment model

syms L K

Taus = 0:1:11;
Ps = 0:0.1:1; %0:0.1:0.9;
doSims = false;
doSims2 =false;
doSims3 = true;

if doSims
Crit = nan(length(Taus),length(Ps));
for t = 1:length(Taus)
    for p = 1:length(Ps)
        
        % find solution for boundary based on Clark.
        % At boundary where |L| = 1, the characteristic equation:
        % L^(Tau+1) - P*L^Tau - f == 0
        % Has a solution e^(i*L)[ or cos(L) + i sin(L) ]
        % f is real, so the imaginary parts of L^(Tau+1) and L^Tau must
        % cancel, and the difference between the real parts must equal f:
        
        Ls = linspace(0,pi/(1+Taus(t)),1e3);
        Ls = Ls(2:end-1); % eliminate the boundary values
        
        S1 = Ps(p)*sin(Taus(t)*Ls);
        S2 = sin((Taus(t)+1)*Ls);
        DS = abs(S1-S2);
        Th = Ls(DS==min(DS));
        
        myfun = @(L) Ps(p)*sin(Taus(t)*L) - sin((Taus(t)+1)*L);
        Th = fzero(@(L) myfun(L), Th);
        
        % Now plug into equation for real part:
        Crit(t,p) = cos((Taus(t)+1)*Th) - Ps(p)*cos(Taus(t)*Th);
         
        
    end % end Ps
end % end Taus


Crit(:,1) = -1; % constrain the values for semelparous univoltine case to be -1
save Clarksims1
else
load Clarksims1
end % end doSims

figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 18])

subplot(3,1,1)
pcolor(Ps,fliplr(Taus),(Crit)); shading flat
colormap(1-gray(length(0:0.2:2)))
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0.5:1:11.5,'yticklabel',10:-1:0)
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'fontsize',10,...
    'ytick',-4:0.2:0)
axis square
set(gca,'fontsize',10)
xlabel(gca,'Adult survival (p)','fontsize',12)
ylabel(gca,'Time to maturity (t, years)','fontsize',12)

% Botsford (1992) version using influence functions:
Taus = 0:1:11;
Ps = 0:0.1:1; %0:0.1:0.9;
As = 1:10;
Kcrit = nan(length(Taus),length(Ps));

if doSims2
for t = 1:length(Taus)
    for p = 1:length(Ps)

        
      %  Now take Botsford influence function approach:
      if Ps(p) > 0
      Infl = (1-Ps(p))*Ps(p).^(As-Taus(t));
      Infl(As<Taus(t)) = 0;
      else
      Infl = 0*As;
      Infl(1)=1;
      end
      

      
      Infl = Infl/sum(Infl);
      %X = sum(Infl.*(1+K).*L.^-As);
      X = sum(Infl.*(K).*L.^-As);
      S = solve(X==1,L,'ReturnConditions',true);
      Ks = linspace(-10,-1,100);
      Ks = Ks(1:end-1);
      
      for k = 1:length(Ks)
      SS_tmp = double(max(abs(vpa(subs(S.L,K,Ks(k))))));
    %  if length(SS_tmp)~=1; keyboard; end
      if isempty(SS_tmp); SS_tmp = NaN; end
      SS(k) = SS_tmp;
      end
      SS = abs(SS-1);
      Kind = find(SS==nanmin(SS),1);
      size(Kind);
      if isempty(Kind)
          Kcrit(t,p) = NaN;
      else
      Kcrit(t,p) = Ks(Kind);
      end
      
               
        
    end % end Ps
end % end Taus
save Clarksims2
else
load Clarksims2
end % end doSims


subplot(3,1,2)
pcolor(Ps,fliplr(Taus),log10(abs(Kcrit))); shading flat; 
colormap((gray(length(0:0.2:2))))
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',0.5:1:11.5,'yticklabel',10:-1:0)
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'fontsize',10)
set(ch,'ytick',log10(2:2:8),'yticklabel',-2:-2:-8)
axis square
set(gca,'fontsize',10)
xlabel(gca,'Adult survival (p)','fontsize',12)
ylabel(gca,'Time to maturity (t, years)','fontsize',12)

doSims3 = false;
if doSims3
    Ks = -1:-1:-10;
    Tau = 2;
    disp(Tau)
for k = 1:length(Ks)
    for p = 1:length(Ps)

        
      %  Now take Botsford influence function approach:
      if Ps(p) > 0
      Infl = (1-Ps(p))*Ps(p).^(As-Tau);
      Infl(As<Tau) = 0;
      else
      Infl = 0*As;
      Infl(1)=1;
      end
          
      Infl = Infl/sum(Infl);
   %   X = sum(Infl.*(1+Ks(k)).*L.^-As);
      X = sum(Infl.*Ks(k).*L.^-As);
      S = solve(X==1,L,'ReturnConditions',true);
     % keyboard
      
     Lam_tmp = max(abs(double((S.L))));
     if isempty(Lam_tmp); Lam_tmp = NaN; end

      Lam(k,p) = Lam_tmp;
      
   
    end % end Ps
end % end Taus
save Clarksims3
else
load Clarksims3
end % end doSims

subplot(3,1,3)
%hold on
pcolor(Ps(2:end),Ks,Lam(:,2:end)); shading flat
%contour(Ps(2:end),Ks,Lam(:,2:end),1);
colormap(1-gray(length(0:0.2:2)))
set(gca,'xtick',0:0.2:1)
set(gca,'ytick',-10:2:0)
ch = colorbar;
set(gca,'tickdir','out','ticklength',[0.015 0.015])
set(ch,'tickdir','out','ticklength',0.015,'fontsize',10)
 %   'ytick',-4:0.2:0)
axis square
set(gca,'fontsize',10)
xlabel(gca,'Adult survival (p)','fontsize',12)
ylabel(gca,'Elasticity (K)','fontsize',12)
      