function bevholt_plaice

% Generate the YPR isopleth figure for North Sea Plaice from B&H 1957

% Dummy data at first
dx = 1;
Ages = 0:dx:15; % age classes
Fs = 0:0.01:2; % values of F to use

% life-history parameters
Linf = 68.5;
k = 0.095;
t0 = -0.815;
M = 0.1;
Amat = 3; % rough estimate from fishbase

Lcs = 4:dx:15; % age-of-entry values to use

Ages = round(Ages,1);
Lcs = round(Lcs,1);

B = Linf*(1 - exp(-k.*(Ages-t0))); % Length-at-age
W = 0.00892*B.^3; % Weight-at-age


for f = 1:length(Fs)
    for a = 1:length(Lcs)
        
% Survival-to-age
     S1 = exp(-M.*Ages(Ages<=Lcs(a))) ;
     S = [S1(end); S1(end).*exp(-(M + Fs(f)).*(Ages(Ages>Lcs(a))'-Lcs(a)))];
     F_frac = Fs(f)./(M + Fs(f)); % proportion of mortality due to fishing
    
     % Yield-at-age (the number that die, times biomass, times proportion
     % fished)
     Y(a,f) = sum(S(:).*(1-exp(-(M+Fs(f)))).*W(Ages>=Lcs(a))'.*F_frac);
     
     % Eggs per recruit
     Surv = [S1(:); S(2:end)];
     Biomass = Surv(:).*W';
     EPR(a,f) = sum(Biomass(Ages>=Amat)); 
     
     % age-specific survival and yield
     SS(a,f,:) = Surv;
     YY{a,f,:} = S(:).*(1-exp(-(M+Fs(f)))).*W(Ages>=Lcs(a))'.*F_frac;
     
        
        if Fs(f)>0
         %   keyboard
        end
        
    end % end Lcs
end % end Fs


% Fig 1: Plot survival & biomass per age, then YPR
% Fig 2: isopleths
figure(1)
clf
set(gcf,'units','cent','position',[10 10 9 11])
sh(1) = subplot(2,1,1);
hold on

As = [2, 7, 11];
FF = Fs==0.4;
for i = 1:3
    SSS = squeeze(SS(As(i),FF,:));
    plot(Ages,SSS,'k-')
end

Pos = get(gca,'position');
axes('position',Pos,'color','none','yaxislocation','right',...
    'tickdir','out','ticklength',[0.02 0.02],'xtick',[])
hold on
plot(Ages,W,'k-')

sh(2) = subplot(2,1,2);
hold on

for i = 1:3
    YYY = squeeze(YY{As(i),FF,:});
    plot(Ages(Ages>=Lcs(As(i))),YYY,'k-')
    sum(YYY)
end
xlim([0 15])


% Fig 2: isopleths
figure(2)
clf
set(gcf,'units','cent','position',[10 10 9 20])


sh(3) = subplot(3,1,1);
hold on
%keyboard
MaxY = max(Y,[],2);
MaxYY = max(Y,[],1);

for a = 1:length(Lcs)
    FmaxY(a) = Fs(Y(a,:)==MaxY(a));
end

for f = 2:length(Fs)
    FmaxYY(f) = Lcs(Y(:,f)==MaxYY(f));
end
contourf(Fs,Lcs,Y,10,'color','k')
plot(FmaxY,Lcs,'ko','markeredgecolor',[1 1 1])
ID = 1:10:200;
plot(Fs(ID),FmaxYY(ID),'kd','markeredgecolor',[1 1 1])
ylim([4 15])
Col = repmat( linspace(0,0.9,10)',[1,3]);
%colormap(1-gray(10))
colormap(1-Col)
ch(1) = colorbar;

sh(4) = subplot(3,1,2);
Profit = Y - 100*repmat(Fs,[length(Lcs),1]);
Profit = max(Profit,0);
contourf(Fs,Lcs,Profit,10,'color','k')
ch(2) = colorbar;

sh(5) = subplot(3,1,3);
EPR = EPR/max(EPR(:));
contourf(Fs,Lcs,EPR,10,'color','k')
ch(3) = colorbar;


set(sh(:),'tickdir','out','ticklength',[0.02 0.02],...
    'xcolor','k','ycolor','k')
set(ch(:),'tickdir','out','ticklength',0.015,...
    'ycolor','k')
set(sh(3:5),'xtick',0:0.4:2, 'ytick',4:2:16)
ylabel(sh(4),'Age of first capture (a_c, y)','fontsize',12)
xlabel(sh(5),'Fishing rate (F, y^-^1)','fontsize',12)

set(sh(1:2),'xtick',0:5:15)
set(sh(1),'ytick',0:0.2:1)
%set(sh(2),'ytick',0:0.1:1)
xlabel(sh(2),'Age (y)','fontsize',12)
ylabel(sh(1),'Survival','fontsize',12)
ylabel(sh(2),'Yield-at-age (g)','fontsize',12)

