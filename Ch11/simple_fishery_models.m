function simple_fishery_models

% Some plots of simple fishery model results

% 1) Shaffer's MSY plot
figure(1)
clf
set(gcf,'units','cent','position',[5 30 9 18])

% Parameters reported by Schaefer for yellowfin tuna:
Binfs = 2.9548e8;
ks = 2.5681;
qs = 3.806e-5;

% Effort (in clipper days)
E= linspace(0,68e3,1e3);


for b = 1:length(Binfs)
for kk = 1:length(ks)
for qq = 1:length(qs)
Binf = Binfs(b);
k = ks(kk);
q = qs(qq);

% Yield as a function of Effort:
Y = Binf.*q.*E - (q.^2.*Binf./k).*E.^2;
Y = max(Y,0);

% Yield per recruit, Bev-Holt method:
F = E.*0.69/25000;
M = 0.8;
tc = 1.5; tr = 1.5; tm = 10; t0= 0.85;
Winf = (0.022.*176.^2.94)./1e3;
Winf = 218;
k = 0.6;
Omega = [1 -3 3 -1];

for n = 1:4
YPR_t(:,n) = F.*exp(-M.*(tc-tr)).*Winf.*Omega(n).*exp(-(n-1).*k.*(tc-t0))...
         .*(1 - exp(-(F + M + (n-1)*k).*(tm-tc)))./(F+M+(n-1)*k);
end
YPR = sum(YPR_t,2);

sh(1)=subplot(2,1,1);
hold on
plot(E,Y,'k')
set(gca,'xlim',[0 68e3],'ylim',[0 2.2e8])
set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k',...
    'ytick',0:0.5e8:2e8)

Pos = get(gca,'position');
Pos(1) = 0.16;
Pos(3) = 0.68;
axes('position',Pos)
plot(E,YPR,'k--')
set(gca,'color','none','yaxislocation','right')

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k')
set(gca,'xlim',[0 68e3],'ylim',[0 22])
ylabel('Yield per recruit, YPR (lbs)','fontsize',12)

set(sh(1),'position',Pos)

axes(sh(1))
xlabel('Fishing effort, E (10^4 days of fishing by class-4 clippers)','fontsize',12)
ylabel('Yield, Y (10^8 lbs)','fontsize',12)
    
% Separate figure: extract stock-recruit relationship.
% Get recruits by dividing Y by Y/R:
R = Y(:)./YPR(:);
% Get stock by dividing Yield/F by q
S = Y(:)./E(:)./qs;


% Plot B: CPUE vs. Yield, to estimate model parameters from data:
subplot(2,1,2)
hold on
set(gca,'position',[Pos(1),0.11,Pos(3:4)])
plot(E,Y./E,'k-')

% some isopleths
Ys = [50 162 250].*1e6;
for y = 1:length(Ys)
    
    plot(E,Ys(y)./E,'k--')
    
end
 

% Load in data
D = importdata('Schaefer_YFT_data.csv');
D = D.data;
plot(D(:,2)*1e3,D(:,3)*1e3,'marker','o','color','k')

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k')
set(gca,'ylim',[0 14e3],'xlim',[0 68e3])


end
end
end

xlabel('Fishing effort, E (10^4 days of fishing by class-4 clippers)','fontsize',12)
ylabel('Catch per unit effort, Y/E (10^3 lbs)','fontsize',12)

% Figure 2: back-calculated stock-recruit curve
figure(2)
clf
set(gcf,'units','cent','position',[20 30 9 9])

plot(S,R,'k-')

xlabel('Stock biomass, B (10^8 lbs)','fontsize',12)
ylabel('Number of recruits (10^6)','fontsize',12)

set(gca,'tickdir','out','ticklength',[0.015 0.015],...
    'xcolor','k','ycolor','k')
set(gca,'ylim',[0 12e6],'xlim',[0 3e8])
axis square

